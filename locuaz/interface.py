from pathlib import Path
from typing import Union, List, Set
import warnings

import freesasa
import MDAnalysis as mda
from Bio.SeqUtils import seq1

from molecules import PDBStructure
from fileutils import FileHandle


def get_interfacing_residues(pdb_input: Union[PDBStructure, FileHandle, Path], chainIDs: List[str]) -> Set[int]:
    """
    get_interfacing_residues(): use freesasa to get the resSeq of the binder residues that are in contact with
    the target. These can be used to guide the choice of the next mutated position.
    Args:
        pdb_input: input PDB.
        chainIDs: only residues belonging to these chains will be reported.

    Returns:
        Set[int]: set of resSeqs from the binder that lie on the interface
    """
    # Remove solvent
    pdb_path = Path(pdb_input)
    u = mda.Universe(str(pdb_path))
    temp_pdb = Path(pdb_path.parent, "temp.pdb")
    u.select_atoms("not resname SOL and not resname WAT").write(str(temp_pdb))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        structs = freesasa.structureArray(str(temp_pdb),
                                      {"separate-chains": False, "chain-groups": ''.join(chainIDs)})

    sasa_whole = freesasa.calc(structs[0])
    sasa_binder = freesasa.calc(structs[1])
    residuos_whole = sasa_whole.residueAreas()
    residuos_binder = sasa_binder.residueAreas()
    interfacing_resis: Set[int] = set()
    for chainID in set(chainIDs):
        for (resnum, sasa_whole), (_, sasa_binder) in zip(residuos_whole[chainID].items(),
                                                          residuos_binder[chainID].items()):
            sasa_diff = sasa_binder.total - sasa_whole.total
            if sasa_diff > 2.:
                interfacing_resis.add(int(resnum))

    # Remove temporaries
    temp_pdb.unlink()

    return interfacing_resis

def get_freesasa_residues(pdb_input: Union[PDBStructure, FileHandle, Path], chainIDs: List[str]) -> Set[int]:
    pdb_path = Path(pdb_input)
    u = mda.Universe(str(pdb_path))
    temp_pdb = Path(pdb_path.parent, "temp.pdb")
    u.select_atoms("not resname SOL and not resname WAT").write(str(temp_pdb))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        structs = freesasa.structureArray(str(temp_pdb),
                                      {"separate-chains": False, "chain-groups": ''.join(set(chainIDs))})
    freesasa_resis = {(int(structs[1].residueNumber(i)), seq1(structs[1].residueName(i))) for i in
                      range(structs[1].nAtoms())}
    temp_pdb.unlink()

    return freesasa_resis