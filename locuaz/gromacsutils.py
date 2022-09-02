from pathlib import Path
from collections.abc import Sequence
from shutil import SameFileError
from typing import Dict
from fileutils import FileHandle, update_header, copy_to
from molecules import (
    AbstractComplex,
    GROComplex,
    PDBStructure,
    get_gro_ziptop_from_pdb,
    copy_mol_to,
    get_tpr,
)
from biobb_md.gromacs.gmxselect import Gmxselect
from biobb_analysis.gromacs.gmx_trjconv_str import gmx_trjconv_str, GMXTrjConvStr
from biobb_md.gromacs.solvate import Solvate
from biobb_md.gromacs.grompp import Grompp
from biobb_md.gromacs.genion import Genion
from biobb_md.gromacs.pdb2gmx import Pdb2gmx
from primitives import launch_biobb
import logging


def write_non_overlapping_ndx(
    pdb_path: Path,
    ndx_fn_in: FileHandle,
    resSeq: int,
    *,
    dist_threshold: float = 0.3,
    gmx_bin: str = "gmx",
) -> FileHandle:

    # First, get the number of waters that overlap
    ndx_fn_wat_out = str(pdb_path.parent / "overlapping.ndx")
    select_overlapping_waters = Gmxselect(
        input_structure_path=str(pdb_path),
        input_ndx_path=str(ndx_fn_in),
        output_ndx_path=str(ndx_fn_wat_out),
        properties={
            "gmx_path": gmx_bin,
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
    ndx_fn_out = str(pdb_path.parent / "non_overlapping.ndx")
    select_not_overlapping = Gmxselect(
        input_structure_path=str(pdb_path),
        input_ndx_path=str(ndx_fn_in),
        output_ndx_path=str(ndx_fn_out),
        properties={
            "gmx_path": gmx_bin,
            "selection": f"not (same residue as resname SOL and "
            f"within {dist_threshold} of (group binder and resid {resSeq})) "
            "and (not name CL) and (not name NA)",
        },
    )
    launch_biobb(select_not_overlapping)

    ndx = FileHandle(ndx_fn_out)
    update_header(ndx, f"[ non_overlapping ]\n")

    wat_count = wat_atoms // 3
    logging.info(
        f"Removed {wat_count} that were within {dist_threshold}nm from "
        "the mutated residue. Will replace them later."
    )

    return ndx, wat_count


def remove_overlapping_waters(
    config: Dict, complex: AbstractComplex, overlapped_resSeq: int
) -> AbstractComplex:
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
    # Get the index file with all the atoms except the waters that overlap
    # with the newly mutated residue
    ndx, wat_count = write_non_overlapping_ndx(
        complex.pdb.file.path,
        complex.ndx,
        overlapped_resSeq,
        gmx_bin=config["md"]["gmx_bin"],
    )

    # Now, write the PDB without the overlapped waters
    nonwat_pdb_fn = Path(complex.dir) / ("nonwat_" + config["main"]["name"] + ".pdb")

    trjconv = GMXTrjConvStr(
        input_structure_path=str(complex.pdb.file.path),
        input_top_path=str(complex.tpr.file.path),
        input_index_path=str(ndx.path),
        output_str_path=str(nonwat_pdb_fn),
        properties={
            "gmx_path": config["md"]["gmx_bin"],
            "selection": "non_overlapping",
        },
    )
    launch_biobb(trjconv)

    # Get .gro and .top for the solvation
    nonwat_pdb, nonwat_gro, nonwat_top = get_gro_ziptop_from_pdb(
        pdb=PDBStructure.from_path(nonwat_pdb_fn),
        target_chains=complex.top.target_chains,
        binder_chains=complex.top.binder_chains,
        water_type=config["md"]["water_type"],
        force_field=config["md"]["force_field"],
        gmx_path=config["md"]["gmx_bin"],
        add_ions=False,
    )

    # pdb2gmx changes the name of the solvent to HOH, for some reason.
    # nonwat_gro.file.replace_text("HOH", "SOL")

    # Re-add waters to keep the N of the system constant
    wet_pdb_fn = Path(complex.dir) / ("wet_" + config["main"]["name"] + ".pdb")
    if wat_count != 0:
        logging.info(f"Adding {wat_count} water molecules.")

        solvatador = Solvate(
            input_solute_gro_path=str(nonwat_pdb),
            input_top_zip_path=str(nonwat_top),
            output_gro_path=str(wet_pdb_fn),
            # actually, the top file is modified inplace by gmx solvent, so this
            # parameter doesn't make much sense.
            output_top_zip_path=str(nonwat_top),
            properties={
                "gmx_path": config["md"]["gmx_bin"],
                "box_type": config["md"]["box_type"],
                "dev": f"-maxsol {wat_count}",
            },
        )
        launch_biobb(solvatador)
    else:
        copy_to(nonwat_pdb.file, complex.dir, wet_pdb_fn.name)
        logging.info(
            f"Not addying any water molecules.{nonwat_pdb} and {wet_pdb_fn} are the same files."
        )

    # Get .gro and .top for ion addition
    wet_pdb, wet_gro, wet_top = get_gro_ziptop_from_pdb(
        pdb=PDBStructure.from_path(wet_pdb_fn),
        target_chains=complex.top.target_chains,
        binder_chains=complex.top.binder_chains,
        water_type=config["md"]["water_type"],
        force_field=config["md"]["force_field"],
        gmx_path=config["md"]["gmx_bin"],
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
            "gmx_path": config["md"]["gmx_bin"],
            "neutral": True,
            "concentration": 0.0,
        },
    )
    launch_biobb(genio)

    # Finally, build the Complex
    new_complex = GROComplex.from_gro_zip(
        name=config["main"]["name"],
        input_dir=complex.dir,
        target_chains=complex.top.target_chains,
        binder_chains=complex.top.binder_chains,
        gmx_bin=config["md"]["gmx_bin"],
    )

    return new_complex
