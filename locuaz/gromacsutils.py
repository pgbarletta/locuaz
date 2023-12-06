import logging
import shutil as sh
import subprocess as sp
import warnings
from pathlib import Path
from typing import Dict, Tuple, Optional, Set

import numpy as np
import MDAnalysis as mda
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO
from biobb_analysis.gromacs.gmx_image import GMXImage
from biobb_analysis.gromacs.gmx_trjconv_str import GMXTrjConvStr
from biobb_gromacs.gromacs.genion import Genion
from biobb_gromacs.gromacs.gmxselect import Gmxselect
from biobb_gromacs.gromacs.solvate import Solvate

from locuaz.complex import AbstractComplex, GROComplex
from locuaz.fileutils import FileHandle, update_header, copy_to
from locuaz.molecules import (
    PDBStructure,
    XtcTrajectory,
    get_tpr,
)
from locuaz.moleculesutils import get_gro_ziptop_from_pdb, fix_wat_naming
from locuaz.primitives import AA_MAP
from locuaz.primitives import launch_biobb


def image_traj(
    cpx: GROComplex, out_trj_fn: Path, use_tleap: bool = False, gmx_bin: str = "gmx"
) -> XtcTrajectory:
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
            "cluster_selection": "Protein",
            "output_selection": "Protein",
            "ur": "compact",
            "pbc": "cluster",
            "center": False,
        },
    )
    launch_biobb(cluster)

    # Use MDAnalysis to get a good reference frame for -pbc nojump
    orig_u = mda.Universe(str(cpx.tpr), str(cpx.tra))
    #### TODO: remove when mdanalysis is updated
    if orig_u.atoms.segids[0][:3] == "seg":
        orig_u.segments.segids = np.array(
            [segid.split("_")[-1] for segid in orig_u.segments.segids]
        )
    #### TODO: remove when mdanalysis is updated
    orig_u.add_TopologyAttr("chainID", orig_u.atoms.segids)

    # First, get a PDB with the same topology as `cluster_trj`.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        orig_pdb = Path(wrk_dir, "orig.pdb")
        # Selection
        orig_u.select_atoms(cpx.top.selection_complex).write(str(orig_pdb))
    u = mda.Universe(str(orig_pdb), str(cluster_trj))
    # noinspection all
    u.trajectory[2]  # type: ignore
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
        },
    )
    launch_biobb(center)

    # Finally, get the last frame of the trajectory and leave it in case the user wants to check it.
    out_pdb_fn = Path(out_trj_fn.parent, out_trj_fn.stem + ".pdb")
    u = mda.Universe(str(orig_pdb), str(out_trj_fn))
    # noinspection all
    u.trajectory[-1]
    if not use_tleap:
        # Staggered resSeq
        for s in u.segments:
            if s.segid in {"", "X"}:
                break
            s.residues.resids = np.array(range(1, len(s.residues) + 1))
    # Continuous resSeq
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        u.atoms.write(str(out_pdb_fn))

    # Remove temporary files
    whole_trj.unlink()
    cluster_trj.unlink()
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


def remove_overlapping_waters(
    complex: GROComplex, config: Dict, overlapped_resSeq: int
) -> GROComplex:
    """DEPRECATED remove_overlapping_waters takes a complex and removes all the waters that
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
        gmx_bin="gmx",
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
            "binary_path": "gmx",
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
                "binary_path": "gmx",
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
    )

    # Build a temporary .tpr for genion
    wet_tpr = get_tpr(gro=wet_gro, top=wet_top, gmx_bin="gmx")

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
            "binary_path": "gmx",
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
        gmx_bin="gmx",
    )

    return new_complex


def fix_gromacs_pdb(
    pdb_in_fn: Path,
    pdb_out_fn: Path,
    *,
    new_chainID: Optional[str] = None,
    allowed_nonstandard_residues: Optional[Set] = None,
) -> Path:
    """
    uses Bio to add TER between chains, END at the end, and manually make the resSeq numbers continuous.
    Since this changes the PBD numbering it should only be used for structures to be used for scoring.
    For ex., Haddock requires chains to have unique resSeq numbers.

    Parameters
    ----------
    pdb_in_fn : Path
        input PDB
    pdb_out_fn : Path
        output PDB
    new_chainID : Optional[str]
        For assigning a new chain ID to all atoms in the output PDB
    allowed_nonstandard_residues : Optional[set]
        Any residue not present in AA_MAP or in this set will be discarded from the output PDB.
        Useful when scoring complexes with ligands.
    Returns
    -------
    the output PDB : Path
        the PDB will be written with Biopython's PDB tools.
    """
    if not allowed_nonstandard_residues:
        allowed_nonstandard_residues = set()
    parsero = PDBParser(QUIET=True)
    pdb = parsero.get_structure("foo", file=pdb_in_fn)

    # Modifying resSeq, chainID and mapping resnames to standard resnames.
    resSeq = 1
    for resi in pdb.get_residues():
        juan, _, roman = resi._id
        try:
            # Rename non-standard residue to its standard equivalent
            resi.resname = AA_MAP[resi.resname]
        except KeyError:
            if resi.resname not in allowed_nonstandard_residues:
                # Not even an amino acid, removing it.
                chain = resi.get_parent()
                chain.detach_child(resi.id)
                continue
        resi._id = (juan, resSeq, roman)
        resSeq += 1
    if new_chainID is not None:
        for chain in pdb.get_chains():
            chain._id = new_chainID

    io = PDBIO()
    io.set_structure(pdb)
    io.save(str(pdb_out_fn))

    return pdb_in_fn


def remove_overlapping_solvent(
    overlapped_pdb: PDBStructure,
    overlapped_resSeq: int,
    nonoverlapped_out_pdb: Path,
    log: logging.Logger,
    *,
    cutoff: float = 3,
    use_tleap=False,
) -> PDBStructure:
    """
    remove_overlapping_solvent() is overlapped with later stages. It may leave a system without ions,
    but not in a state that would imply discontinuity between ions and solvent, when ions are added
    on a later step.
    Thus, it removes all ions, relying on later steps (the building of Complex) to add them if necessary.
    """
    pdb_in_fn = Path(overlapped_pdb)
    u = mda.Universe(str(pdb_in_fn))

    overlapped_wat_atoms = u.select_atoms(
        f"(around {cutoff} resnum {overlapped_resSeq}) and (resname WAT or resname SOL)"
    ).residues.atoms
    nwats_atm = len(overlapped_wat_atoms)
    assert (
        nwats_atm % 3
    ) == 0, f"Invalid number of overlapped water atoms: {nwats_atm}"
    nwats = len(overlapped_wat_atoms) // 3

    ions = u.select_atoms(
        "name CL or name Cl or name NA or name Na or type CL or type Cl or type NA or type Na"
    )
    nions = len(ions)

    nonoverlapped_nonwat = "init_nonoverlapped_nonwat.pdb"
    nonoverlapped_nonwat_fn = Path(pdb_in_fn.parent, nonoverlapped_nonwat)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        (u.atoms - overlapped_wat_atoms - ions).write(str(nonoverlapped_nonwat_fn))  # type: ignore

    # insert waters
    nonoverlapped = "init_nonoverlapped.pdb"
    nonoverlapped_fn = Path(pdb_in_fn.parent, nonoverlapped)
    if (nwats + nions) == 0:
        log.info(
            f"Not addying any water molecules.{nonoverlapped_nonwat} and "
            f"{nonoverlapped} are the same files."
        )
        sh.copy(nonoverlapped_nonwat_fn, nonoverlapped_fn)
        fix_wat_naming(nonoverlapped_fn, nonoverlapped_out_pdb, use_tleap=use_tleap)
    else:
        log.info(
            f"Removed {nwats} water molecules within {cutoff}A from the mutated residue. "
            "Will replace them later."
        )
        comando_solv = f"gmx solvate -cp {nonoverlapped_nonwat} -o {nonoverlapped} -maxsol {nwats + nions}"
        log.info(
            f"Adding {nwats} water molecules, and another {nions} that will be replaced by ions later."
        )
        p = sp.run(
            comando_solv,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            cwd=nonoverlapped_fn.parent,
            shell=True,
            text=True,
        )
        try:
            fix_wat_naming(nonoverlapped_fn, nonoverlapped_out_pdb, use_tleap=use_tleap)
        except Exception as e:
            raise Exception(
                f"gmx solvate failed when inserting waters: {p.stdout} \n {p.stderr}"
            ) from e

    return PDBStructure.from_path(nonoverlapped_fn)
