from pathlib import Path
from typing import Tuple, Optional, Union

import MDAnalysis as mda
from biobb_analysis.gromacs.gmx_trjconv_str import GMXTrjConvStr

from complex import GROComplex
from molecules import PDBStructure
from primitives import launch_biobb


def split_solute_and_solvent_old(
        cpx: GROComplex, gmx_bin: str
) -> Tuple[PDBStructure, PDBStructure]:
    """prepare_old_iter extract 2 PDBs from an input pdb, one with the protein
    and the other with the water and ions.

    Args:
        cpx (Complex): a complex object with a PDB and a TPR file.
        gmx_bin:

    Returns:
        Tuple[PDBStructure, PDBStructure]: solute and solvent+ions
    """

    # Protein
    nonwat_pdb_fn = Path(cpx.dir) / ("nonwat_" + cpx.name + ".pdb")
    get_protein = GMXTrjConvStr(
        input_structure_path=str(cpx.pdb.file.path),
        input_top_path=str(cpx.tpr.file.path),
        input_index_path=str(cpx.ndx.path),
        output_str_path=str(nonwat_pdb_fn),
        properties={"binary_path": gmx_bin, "selection": "Protein"},
    )
    launch_biobb(get_protein)
    nonwat_pdb = PDBStructure.from_path(nonwat_pdb_fn)

    # Water and ions
    wation_pdb_fn = Path(cpx.dir) / ("wation_" + cpx.name + ".pdb")
    get_water_ions = GMXTrjConvStr(
        input_structure_path=str(cpx.pdb.file.path),
        input_top_path=str(cpx.tpr.file.path),
        input_index_path=str(cpx.ndx.path),
        output_str_path=str(wation_pdb_fn),
        properties={"binary_path": gmx_bin, "selection": "Non-Protein"},
    )
    launch_biobb(get_water_ions)

    wation_pdb = PDBStructure.from_path(wation_pdb_fn)

    return nonwat_pdb, wation_pdb


def split_solute_solvent(
    pdb_in: Union[PDBStructure, Path],
    *,
    selection_complex: Optional[str] = None,
    selection_wations: Optional[str] = None,
) -> Tuple[PDBStructure, PDBStructure]:

    pdb_in_path = Path(pdb_in)
    u = mda.Universe(str(pdb_in_path))

    # Protein
    nonwat_pdb_fn = Path(pdb_in_path.parent, "init_nonwat.pdb")
    if selection_complex:
        u.atoms.select_atoms(selection_complex).write(str(nonwat_pdb_fn))  # type: ignore
    else:
        u.atoms.select_atoms("protein").write(str(nonwat_pdb_fn))  # type: ignore

    # Water and ions
    wation_pdb_fn = Path(pdb_in_path.parent, "init_wation.pdb")
    if selection_wations:
        u.atoms.select_atoms(selection_wations).write(str(wation_pdb_fn))  # type: ignore
    else:
        u.atoms.select_atoms("not protein").write(str(wation_pdb_fn))  # type: ignore

    return PDBStructure.from_path(nonwat_pdb_fn), PDBStructure.from_path(wation_pdb_fn)
