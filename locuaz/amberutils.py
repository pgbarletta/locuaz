from pathlib import Path
from collections import deque
from collections.abc import Iterable
from typing import Tuple, Union, Optional, List, Generator, Deque, Set
import subprocess as sp
from zipfile import ZipFile
import shutil as sh
from io import TextIOWrapper
import warnings
import sys

import parmed as pmd
from pdb4amber import AmberPDBFixer
import MDAnalysis as mda

from fileutils import FileHandle, DirHandle
from molecules import PDBStructure, GROStructure, ZipTopology
from primitives import ext

warnings.filterwarnings("ignore")


def chunk(f: TextIOWrapper) -> Generator:
    linea = " "
    while linea != "":
        text: List[str] = []
        while True:
            cursor = f.tell()
            linea = f.readline()
            mol_start = linea[:16] == "[ moleculetype ]"
            last_mol_end = linea[:10] == "[ system ]"
            if (mol_start or last_mol_end) and len(text) > 0:
                f.seek(cursor)
                break
            elif linea == "":
                break
            else:
                text.append(linea)
        yield text


def create_tleap_script(local_dir: Union[DirHandle, Path], name: str) -> FileHandle:
    in_tleap = Path(local_dir, "tleap")
    out_contents = []

    with open(in_tleap, "r") as f:
        mol: Optional[str] = None
        for linea in f:
            if "loadpdb" in linea:
                mol = linea.split("=")[0].strip()
                out_contents.append(f"{mol} = loadpdb {ext(name, 'pdb')}")
            elif linea[0:13] == "saveamberparm":
                assert mol, "Bad tleap script. Please follow the expected format."
                out_contents.append(
                    f"saveamberparm {mol} pre_{ext(name, 'prmtop')} pre_{ext(name, 'rst7')}"
                    "\n"
                )
            else:
                out_contents.append(linea)

    out_tleap = Path(local_dir, f"tleap_{name}")
    with open(out_tleap, "w") as f:
        for line in out_contents:
            f.write(line)

    # Remove temporary tleap script.
    in_tleap.unlink()

    return FileHandle(out_tleap)


def fix_pdb(in_pdb: Union[PDBStructure, Path], out_pdb_path: Path) -> PDBStructure:

    pdb_path = Path(in_pdb)

    pdbfixer = AmberPDBFixer(str(pdb_path))
    sslist, _ = pdbfixer.find_disulfide()
    pdbfixer.rename_cys_to_cyx(sslist)
    pdbfixer.assign_histidine()

    pdbfixer.parm.save(str(out_pdb_path), conect=False, overwrite=True)

    return PDBStructure(FileHandle(out_pdb_path))


def run_tleap(
    tleap_script: Union[FileHandle, Path], name: str
) -> Tuple[PDBStructure, FileHandle, FileHandle]:
    """run_tleap(): runs tleap script and then uses a local PDB to add chainID info
    to the .prmtop and box info to the .rst7. Then it uses both to overwrite the
    PDB to keep everything in sync (in case, for example, that the tleap script
    added ions)


    Args:
        tleap_script (Union[FileHandle, Path]): script that loads a local PDB and
        outputs a .prmtop and a .rst7 file
        name (str): name of the system

    Returns:
        Tuple[PDBStructure, FileHandle, FileHandle]: PDB, .prmtop and .rst7
    """

    local_dir = Path(tleap_script).parent

    p = sp.run(
        f"tleap -f {Path(tleap_script)}",
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        cwd=local_dir,
        shell=True,
        text=True,
    )

    pre_top_path = Path(local_dir, f"pre_{name}.prmtop")
    pre_rst_path = Path(local_dir, f"pre_{name}.rst7")
    pdb_path = Path(local_dir, f"{name}.pdb")
    assert (
        pre_top_path.is_file() and pre_rst_path.is_file() and pdb_path.is_file()
    ), f"---- Possible tleap error: ----\n{p.stdout}"

    amb = pmd.load_file(str(pre_top_path), str(pre_rst_path))
    u: mda.Universe = mda.Universe(str(pdb_path))

    mda_res = len(u.residues)  # type: ignore
    pmd_res = len(amb.residues)
    assert (
        mda_res == pmd_res
    ), f"Mismatch between parmed residues ({pmd_res}) and MDAnalysis residues ({mda_res})."

    # Set box size:
    amb.box = u.dimensions

    # Add chainID info to prmtop
    chainids: List[str] = []
    for mda_res, pmd_res in zip(u.residues, amb.residues):  # type: ignore
        segid = mda_res.segid if mda_res.segid != "" else "X"  # type: ignore
        pmd_res.chain = segid  # type: ignore
        chainids.append(segid)  # type: ignore

    amb.add_flag(
        "RESIDUE_CHAINID",
        "20a4",
        data=chainids,
        comments=["Residue chain ID (chainId) read from PDB " "file; DIMENSION(NRES)"],
    )

    top_path = Path(local_dir, f"{name}.prmtop")
    rst_path = Path(local_dir, f"{name}.rst7")
    amb.save(str(pdb_path), conect=False, overwrite=True)
    amb.save(str(top_path))
    amb.save(str(rst_path))

    pdb = PDBStructure.from_path(pdb_path)
    prmtop = FileHandle(top_path)
    rst = FileHandle(rst_path)

    return pdb, prmtop, rst


def amb_to_gmx(
    name: str,
    pdb: PDBStructure,
    prmtop: Union[FileHandle, Path],
    rst: Union[FileHandle, Path],
    *,
    target_chains: Iterable[str],
    binder_chains: Iterable[str],
    zip_top=True,
) -> Tuple[GROStructure, ZipTopology]:
    local_dir = Path(prmtop).parent
    amb = pmd.load_file(str(prmtop), str(rst))
    amb_pdb = pmd.load_file(str(pdb))
    amb.box = amb_pdb.box

    gro_path = Path(local_dir, f"{name}.gro")
    top_path = Path(local_dir, f"{name}.top")
    amb.save(str(gro_path))
    amb.save(str(top_path))

    gro = GROStructure(FileHandle(gro_path))

    ziptop_handle = fixup_top(
        top_path,
        name,
        target_chains=list(target_chains),
        binder_chains=list(binder_chains),
    )
    ziptop = ZipTopology.from_path_with_chains(
        ziptop_handle.path, target_chains=target_chains, binder_chains=binder_chains
    )
    # ziptop.target_chains = tuple(target_chains)
    # ziptop.binder_chains = tuple(binder_chains)

    return gro, ziptop


def fixup_top(
    top_path: Path,
    name: str,
    *,
    target_chains: Iterable[str],
    binder_chains: Iterable[str],
) -> FileHandle:
    local_dir = top_path.parent
    # Backup topology
    norestr_name = "norestr_" + top_path.name.split(".")[0] + ".top"
    sh.copy(top_path, local_dir / norestr_name)

    # Get the topology of each molecule
    with open(top_path, "r") as file:
        chunker = chunk(file)
        header_top: List[str] = next(chunker)

        mols_top_text: List[List[str]] = []
        for text in chunker:
            mols_top_text.append(text)

        bottom_top: List[str] = mols_top_text.pop()

    # Write out the topology of each molecule in its dedicated itp
    mol_itp_files = {}

    chainIDs: Deque[str] = deque(target_chains + binder_chains)  # type: ignore
    for mol_text in mols_top_text:
        # Assuming name is 'systemN' where N is the chain number.
        # If they change this, my code breaks.
        mol_name_old = mol_text[2][:10].strip()
        if mol_name_old[:6] == "system":
            mol = chainIDs.popleft()
            # mol(chainID) should be just a char, but just in case:
            mol_text[2] = f"{mol}{3: >{21 - len(mol)}}" "\n"
        else:
            mol = mol_name_old

        mol_itp_files[mol] = Path(local_dir, f"{mol}.itp")
        with open(mol_itp_files[mol], "w") as file:
            for line in mol_text:
                file.write(line)

    # Open them with MDA, select the heavy atoms and write the posres_.itp restraints files
    header = """;automatically generated by LOCUAZ

    [ position_restraints ]
    ; atom  type      fx      fy      fz
    """
    wat_header = """automatically generated by LOCUAZ
    [ position_restraints ]
    ;  i funct       fcx        fcy        fcz
    1    1       1000       1000       1000
    """
    ion_text = """     1     1  1000  1000  1000
    """
    mol_posres_files = {}
    protein: Set[str] = set(target_chains + binder_chains)  # type: ignore
    solvent: Set[str] = {"SOL", "WAT"}
    ions: Set[str] = {"NA", "Na", "CL", "Cl"}
    for mol, mol_path in mol_itp_files.items():
        mol_posres_files[mol] = Path(mol_path.parent, f"posre_{mol}.itp")
        assert not mol_posres_files[mol].is_file(), "This shouldn't happen"
        file = open(mol_posres_files[mol], "w")

        if mol in protein:
            file.write(header)

            u = mda.Universe(str(mol_path))
            at = u.select_atoms("not type H and not type H?")
            for a in at:
                file.write(f"{a.id:>6d}     1  1000  1000  1000" "\n")
        elif mol in solvent:
            file.write(wat_header)
        elif mol in ions:
            file.write(header)
            file.write(ion_text)
        else:
            warnings.warn(
                f"Writing {mol}'s topology as if it was an ion. Check if it's ok."
            )
            file.write(header)
            file.write(ion_text)
        file.close()

    # Add posres_.itp include on each molecules.itp
    for mol, mol_path in mol_itp_files.items():
        posres_include_text = f"""; Include Position restraint file
#ifdef POSRES
#include "posre_{mol}.itp"
#endif

"""
        with open(mol_path, "a") as file:
            file.write(posres_include_text)

    includes_top = ["; Include chain topologies" "\n"]
    for mol_path in mol_itp_files.values():
        includes_top.append(f'#include "{mol_path.name}"' "\n")

    # create INCLUDES FOR THE molecules.itp, merge [ header_top ; includes_top ; bottom_top ] and write main top file
    chainIDs: Deque[str] = deque(target_chains + binder_chains)  # type: ignore
    with open(top_path, "w") as file:
        for line in header_top + ["\n"] + includes_top:
            file.write(line)
        for line in bottom_top:
            # rewrite 'system1', 'system2', ... with the proteins chainIDs
            if line[:6] == "system":
                chainID = chainIDs.popleft()
                # chainID shoul be just a char, but just in case:
                line = f"{chainID}{1: >{22 - len(chainID)}}" "\n"
            file.write(line)

    # compress everything on a zip file
    ziptop_path = Path(local_dir, f"{name}.zip")
    with ZipFile(ziptop_path, mode="w") as zf:
        zf.write(top_path, arcname=f"{name}.top")
        for mol, mol_path in mol_itp_files.items():
            zf.write(mol_path, arcname=f"{mol}.itp")
        for mol, posres_path in mol_posres_files.items():
            zf.write(posres_path, arcname=f"posre_{mol}.itp")

    return FileHandle(ziptop_path)


# DEPRECATED
def run_acpype(
    name: str,
    prmtop: Union[FileHandle, Path],
    rst: Union[FileHandle, Path],
    *,
    target_chains: Iterable,
    binder_chains: Iterable,
    zip_top=True,
) -> Tuple[GROStructure, ZipTopology]:
    local_dir = Path(prmtop).parent
    p = sp.run(
        f"acpype -p {prmtop} -x {rst}",
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        cwd=local_dir,
        shell=True,
        text=True,
    )

    error_msg = "Execution of acpype failed."
    assert p.stdout, error_msg

    acpype_dir = Path(local_dir, f"{name}.amb2gmx")
    acpype_gro_path = Path(acpype_dir, f"{name}_GMX.gro")
    gro_path = Path(acpype_dir, f"{name}.gro")
    acpype_top_path = Path(acpype_dir, f"{name}_GMX.top")
    top_path = Path(acpype_dir, f"{name}.top")
    ziptop_path = Path(local_dir, f"{name}.zip")
    try:
        sh.move(acpype_gro_path, gro_path)
        sh.move(acpype_top_path, top_path)

        gro = GROStructure(FileHandle(gro_path))
        with ZipFile(ziptop_path, mode="w") as zf:
            zf.write(top_path, arcname=f"{name}.top")
        ziptop = ZipTopology(FileHandle(ziptop_path))
        ziptop.target_chains = tuple(target_chains)
        ziptop.binder_chains = tuple(binder_chains)

    except Exception as e:
        raise FileNotFoundError(error_msg) from e

    return gro, ziptop