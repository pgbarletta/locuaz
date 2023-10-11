import warnings
from pathlib import Path
from typing import Union, List, Set, Optional, Tuple
import py

import freesasa
import MDAnalysis as mda
from Bio.SeqUtils import seq1

from locuaz.molecules import PDBStructure
from locuaz.fileutils import FileHandle


def get_interfacing_residues(pdb_input: Union[PDBStructure, FileHandle, Path],
                             chainIDs: List[str],
                             probe_radius: float = 1.4,
                             amber_numbering: bool = False) -> Set[int]:
    """
    get_interfacing_residues(): use freesasa to get the resSeq of the binder residues that are in contact with
    the target. These can be used to guide the 2choice of the next mutated position.
    Args:
        pdb_input (Union[PDBStructure, FileHandle, Path]): input PDB.
        chainIDs (List[str]): only residues belonging to these chains will be reported.
        probe_radius: probe_radius for freesasa
        amber_numbering (bool): if True, it will renumber the resSeqs of the input PDB (on another temporary PDB)
        as per Amber numbering scheme, that is, continuous resSeq numbers as opposed to GROMACS numbering
        which starts at 1 on each chain.

    Returns:
        Set[int]: set of resSeqs from the binder that lie on the interface
    """
    # Remove solvent
    pdb_path = Path(pdb_input)
    u = mda.Universe(str(pdb_path))
    temp_pdb = Path(pdb_path.parent, "temp.pdb")
    complex = u.select_atoms("not (resname SOL or resname WAT or resname CL or resname NA or resname Cl or resname Na)")
    # Renumber resSeq when using Amber numbering since GROMACS will have renumbered the PDB
    # to start at 1 on each chain
    if amber_numbering:
        complex.residues.resids = range(1, len(complex.residues.resnums) + 1)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        complex.write(str(temp_pdb))

    # Silence warnings from freesasa
    capture_warnings = py.io.StdCaptureFD(out=True, in_=False)
    structs = freesasa.structureArray(str(temp_pdb),
                                      {"separate-chains": False,
                                       "hetatm": True,
                                       "chain-groups": ''.join(set(chainIDs))})
    capture_warnings.reset()

    parameters = freesasa.Parameters({
        'algorithm': "LeeRichards",
        'probe-radius': probe_radius,
        'n-points': 200,
        'n-slices': 20,
        'n-threads': 1})
    sasa_whole = freesasa.calc(structs[0], parameters)
    sasa_binder = freesasa.calc(structs[1], parameters)

    residuos_whole = sasa_whole.residueAreas()
    residuos_binder = sasa_binder.residueAreas()
    interfacing_resis: Set[int] = set()
    for chainID in set(chainIDs):
        for (resnum, sasa_whole), (_, sasa_binder) in zip(residuos_whole[chainID].items(),
                                                          residuos_binder[chainID].items()):
            sasa_diff = sasa_binder.total - sasa_whole.total
            if sasa_diff > 0.5:
                interfacing_resis.add(int(resnum))

    # Remove temporaries
    temp_pdb.unlink()

    return interfacing_resis


def get_freesasa_residues(pdb_input: Union[PDBStructure, FileHandle, Path], chainIDs: List[str]) -> Set[int]:
    pdb_path = Path(pdb_input)
    u = mda.Universe(str(pdb_path))
    temp_pdb = Path(pdb_path.parent, "temp.pdb")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        u.select_atoms("not resname SOL and not resname WAT").write(str(temp_pdb))

    # Silence warnings from freesasa
    capture_warnings = py.io.StdCaptureFD(out=True, in_=False)
    structs = freesasa.structureArray(str(temp_pdb),
                                      {"separate-chains": False,
                                       "hetatm": False,
                                       "chain-groups": ''.join(set(chainIDs))})
    capture_warnings.reset()

    freesasa_resis = {(int(structs[1].residueNumber(i)), seq1(structs[1].residueName(i))) for i in
                      range(structs[1].nAtoms())}
    temp_pdb.unlink()

    return freesasa_resis


def get_interface_surface(pdb_path: Path, i: Optional[int] = None) -> Tuple[int, Optional[float]]:
    """
    get_interface_surface(): `i` could be the index of the input PBD, in case the function is called
    on a concurrent manner and order is needed when populating an array
    Args:
        pdb_path: input PDB
        i: optional index of the input PBD, it'll be returned as is. Useful in case the function is called on a concurrent manner and order is needed when populating an array

    Returns:
        Tuple[int, float] = i, surface in A^2
    """

    # Silence warnings from freesasa
    capture_warnings = py.io.StdCaptureFD(out=True, in_=False)
    # Assume A=target, B=binder convention
    structs = freesasa.structureArray(str(pdb_path),
                                      {"separate-chains": False,
                                       "hetatm": False,
                                       "chain-groups": "A+B"})
    capture_warnings.reset()

    parameters = freesasa.Parameters({
        'algorithm': "LeeRichards",
        'probe-radius': 1.4,
        'n-points': 200,
        'n-slices': 20,
        'n-threads': 1})
    sasa_whole = freesasa.calc(structs[0], parameters)
    sasa_a = freesasa.calc(structs[1], parameters)
    sasa_b = freesasa.calc(structs[2], parameters)

    return i, (sasa_a.totalArea() + sasa_b.totalArea() - sasa_whole.totalArea())
