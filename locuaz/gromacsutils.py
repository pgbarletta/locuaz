from functools import singledispatch
from pathlib import Path
from typing import Dict, Tuple
import logging
import warnings

import MDAnalysis as mda

from fileutils import FileHandle, update_header, copy_to
from molecules import (
    PDBStructure,
    Trajectory,
    XtcTrajectory,
    get_gro_ziptop_from_pdb,
    get_tpr,
)
from complex import AbstractComplex, GROComplex
from biobb_gromacs.gromacs.gmxselect import Gmxselect
from biobb_analysis.gromacs.gmx_trjconv_str import GMXTrjConvStr
from biobb_gromacs.gromacs.solvate import Solvate
from biobb_gromacs.gromacs.genion import Genion
from biobb_analysis.gromacs.gmx_image import GMXImage
from primitives import launch_biobb


@singledispatch
def image_traj(cpx: AbstractComplex, out_trj_fn: Path, gmx_bin: str) -> Trajectory:
    raise NotImplementedError


@image_traj.register
def _(cpx: GROComplex, out_trj_fn: Path, gmx_bin: str) -> XtcTrajectory:

    wrk_dir = out_trj_fn.parent

    # I have to specify `fit_selection` and `center_selection` even though I'm not
    # fitting or centering anything or else biobb will throw an error.
    whole_trj = Path(wrk_dir, "whole.xtc")
    make_whole = GMXImage(
        input_traj_path=str(cpx.tra),
        input_index_path=str(cpx.ndx),
        input_top_path=str(cpx.tpr),
        output_traj_path=str(whole_trj),
        properties={
            "binary_path": gmx_bin,
            "fit_selection": "Protein",
            "center_selection": "Protein",
            "output_selection": "Protein",
            "ur": "compact",
            "pbc": "whole",
            "center": False,
        },
    )
    launch_biobb(make_whole)

    cluster_trj = Path(wrk_dir, "clustered.xtc")
    cluster = GMXImage(
        input_traj_path=str(whole_trj),
        input_index_path=str(cpx.ndx),
        input_top_path=str(cpx.tpr),
        output_traj_path=str(cluster_trj),
        properties={
            "binary_path": gmx_bin,
            "fit_selection": "Protein",
            "center_selection": "Protein",
            "output_selection": "Protein",
            "ur": "compact",
            "pbc": "cluster",
            "center": False,
        },
    )
    launch_biobb(cluster)

    # Use MDAnalysis to get a good reference frame for -pbc nojump
    orig_u = mda.Universe(str(cpx.tpr), str(cpx.tra))
    # First, get a PDB with the same topology as `cluster_trj`.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        orig_pdb = Path(wrk_dir, "orig.pdb")
        # Selection 'complex' was made using MDA and selecting 'protein'.
        # If that changes in the future, then this'll break. Sorry.
        orig_u.select_atoms("protein").write(str(orig_pdb))
    u = mda.Universe(str(orig_pdb), str(cluster_trj))
    u.trajectory[2]
    cluster_gro = Path(wrk_dir, "clustered.gro")
    u.atoms.write(str(cluster_gro))  # type: ignore

    nojump_trj = Path(wrk_dir, "nojump.xtc")
    fix_jump = GMXImage(
        input_traj_path=str(cluster_trj),
        input_index_path=str(cpx.ndx),
        input_top_path=str(cluster_gro),
        output_traj_path=str(nojump_trj),
        properties={
            "binary_path": gmx_bin,
            "fit_selection": "Protein",
            "center_selection": "Protein",
            "output_selection": "Protein",
            "ur": "compact",
            "pbc": "nojump",
            "center": False,
        },
    )
    launch_biobb(fix_jump)

    center = GMXImage(
        input_traj_path=str(nojump_trj),
        input_index_path=str(cpx.ndx),
        input_top_path=str(cluster_gro),
        output_traj_path=str(out_trj_fn),
        properties={
            "binary_path": gmx_bin,
            "fit_selection": "Protein",
            "center_selection": "Protein",
            "output_selection": "Protein",
            "ur": "compact",
            "center": True,
            "center_selection": "Protein",
        },
    )
    launch_biobb(center)

    # Remove temporary files
    whole_trj.unlink()
    cluster_trj.unlink()
    # orig_zip.unlink()
    orig_pdb.unlink()
    cluster_gro.unlink()
    nojump_trj.unlink()

    return XtcTrajectory(FileHandle(out_trj_fn))


def write_non_overlapping_ndx(
    pdb_path: Path,
    ndx_fn_in: FileHandle,
    resSeq: int,
    *,
    dist_threshold: float = 0.3,
    gmx_bin: str = "gmx",
) -> Tuple[FileHandle, int]:

    # First, get the number of waters that overlap
    ndx_fn_wat_out = str(pdb_path.parent / "overlapping.ndx")
    select_overlapping_waters = Gmxselect(
        input_structure_path=str(pdb_path),
        input_ndx_path=str(ndx_fn_in),
        output_ndx_path=str(ndx_fn_wat_out),
        properties={
            "binary_path": gmx_bin,
            "selection": f"(same residue as resname SOL and "
            f"within {dist_threshold} of (group binder and resid {resSeq}))",
        },
    )
    launch_biobb(select_overlapping_waters)

    wat_atoms = 0
    with open(ndx_fn_wat_out, "r") as file:
        # Discard header
        next(file)
        for linea in file:
            # Count how many words (atom serial numbers) are in each line and
            # substract one because the new line character will be counted as well.
            wat_atoms += len(linea.split(" ")) - 1

    if (wat_atoms % 3) != 0:
        raise RuntimeError(
            "Number of overlapping water atoms is not divisible by 3. "
            "This should never happen."
        )

    # Now, get the negation of the previous selection.
    # And also discard ions (Na and Cl). These will be readded later.
    # This is done to prevent the SOL group from being discontinuous (as water
    # molecules will be readded later), which breaks gmx genion.
    ndx_fn_out = Path(pdb_path.parent, "non_overlapping.ndx")
    select_not_overlapping = Gmxselect(
        input_structure_path=str(pdb_path),
        input_ndx_path=str(ndx_fn_in),
        output_ndx_path=str(ndx_fn_out),
        properties={
            "binary_path": gmx_bin,
            "selection": f"not (same residue as resname SOL and "
            f"within {dist_threshold} of (group binder and resid {resSeq})) "
            "and (not name CL) and (not name NA)",
        },
    )
    launch_biobb(select_not_overlapping)

    ndx = FileHandle(ndx_fn_out)
    update_header(ndx, f"[ non_overlapping ]\n")

    wat_count = wat_atoms // 3

    return ndx, wat_count


@singledispatch
def remove_overlapping_waters(
    complex: AbstractComplex, config: Dict, overlapped_resSeq: int
) -> AbstractComplex:
    raise NotImplementedError


@remove_overlapping_waters.register
def _(complex: GROComplex, config: Dict, overlapped_resSeq: int) -> GROComplex:
    """remove_overlapping_waters takes a complex and removes all the waters that
    are within a certain cutoff from `overlapped_resSeq`. Then it replaces them.
    It does a lot of gymnastics to deal with Gromacs weirdness.

    Args:
        config (Dict): work project's config file
        complex (AbstractComplex): a Gromacs complex
        overlapped_resSeq (int): assumes that this residue belongs to 'binder'

    Returns:
        AbstractComplex: fully loaded complex with no overlaps between water
        and `overlapped_resSeq`
    """
    log = logging.getLogger(config["main"]["name"])
    # Get the index file with all the atoms except the waters that overlap
    # with the newly mutated residue
    ndx, wat_count = write_non_overlapping_ndx(
        complex.pdb.file.path,
        complex.ndx,
        overlapped_resSeq,
        gmx_bin=config["md"]["gmx_bin"],
    )

    log.info(
        f"Removed {wat_count} water molecules within 0.3nm from the mutated residue. "
        "Will replace them later."
    )

    # Now, write the PDB without the overlapped waters
    nonwat_pdb_fn = Path(complex.dir) / ("nonwat_" + config["main"]["name"] + ".pdb")

    trjconv = GMXTrjConvStr(
        input_structure_path=str(complex.pdb.file.path),
        input_top_path=str(complex.tpr.file.path),
        input_index_path=str(ndx.path),
        output_str_path=str(nonwat_pdb_fn),
        properties={
            "binary_path": config["md"]["gmx_bin"],
            "selection": "non_overlapping",
        },
    )
    launch_biobb(trjconv)

    # Get .gro and .top for the solvation
    nonwat_pdb, nonwat_gro, nonwat_top = get_gro_ziptop_from_pdb(
        pdb=PDBStructure.from_path(nonwat_pdb_fn),
        target_chains=complex.top.target_chains,
        binder_chains=complex.top.binder_chains,
        md_config=config["md"],
        add_ions=False,
    )

    # Re-add waters to keep the N of the system constant
    wet_pdb_fn = Path(complex.dir) / ("wet_" + config["main"]["name"] + ".pdb")
    if wat_count != 0:
        log.info(f"Adding {wat_count} water molecules.")

        solvatador = Solvate(
            input_solute_gro_path=str(nonwat_pdb),
            input_top_zip_path=str(nonwat_top),
            output_gro_path=str(wet_pdb_fn),
            # actually, the top file is modified inplace by gmx solvent, so this
            # parameter doesn't make much sense.
            output_top_zip_path=str(nonwat_top),
            properties={
                "binary_path": config["md"]["gmx_bin"],
                "dev": f"-maxsol {wat_count}",
            },
        )
        launch_biobb(solvatador)
    else:
        copy_to(nonwat_pdb.file, complex.dir, wet_pdb_fn.name)
        log.info(
            f"Not addying any water molecules.{nonwat_pdb} and {wet_pdb_fn} are the same files."
        )

    # Get .gro and .top for ion addition
    wet_pdb, wet_gro, wet_top = get_gro_ziptop_from_pdb(
        pdb=PDBStructure.from_path(wet_pdb_fn),
        target_chains=complex.top.target_chains,
        binder_chains=complex.top.binder_chains,
        md_config=config["md"],
        add_ions=False,
    )

    # Build a temporary .tpr for genion
    wet_tpr = get_tpr(gro=wet_gro, top=wet_top, gmx_bin=config["md"]["gmx_bin"])

    # Re-add ions as necessary. SOL group will be continous, so gmx genion
    # won't complain.
    gro_fn = Path(complex.dir) / (config["main"]["name"] + ".gro")
    top_fn = Path(complex.dir) / (config["main"]["name"] + ".zip")
    genio = Genion(
        input_tpr_path=str(wet_tpr),
        input_top_zip_path=str(wet_top),
        output_gro_path=str(gro_fn),
        output_top_zip_path=str(top_fn),
        properties={
            "binary_path": config["md"]["gmx_bin"],
            "neutral": True,
            "concentration": 0.0,
        },
    )
    launch_biobb(genio)

    # Finally, build the Complex
    new_complex: GROComplex = GROComplex.from_gro_zip(
        name=config["main"]["name"],
        input_dir=Path(complex.dir),
        target_chains=complex.top.target_chains,
        binder_chains=complex.top.binder_chains,
        gmx_bin=config["md"]["gmx_bin"],
    )

    return new_complex
