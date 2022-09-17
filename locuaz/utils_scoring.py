import re
from typing import Iterable
import zipfile
from pathlib import Path
import glob

from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO

from fileutils import FileHandle, catenate_pdbs


def extract_pdbs(
    zip_file: Path,
    out_prefix: str,
    *,
    new_chainID: str = None,
    delete_zip_file: bool = True,
):
    """
    zip_file: zip file with a bunch of PDBs named "output{i}.pdb" generated by trjconv
    Extracts the zip file, fixes the PDBs and, optionally, deletes the original zipfile
    """
    zip_file = Path(zip_file)
    current_dir = zip_file.parent
    zipped_pdbs = zipfile.ZipFile(zip_file)
    regex_get_number = r"\d+"
    nframes = 1 + max(
        [
            int(re.findall(regex_get_number, nombre.filename)[0])
            for nombre in zipped_pdbs.filelist
        ]
    )

    with zipped_pdbs as sipesipe:
        sipesipe.extractall(current_dir)

    for i in range(nframes):
        old_pdb = Path(current_dir / (f"output{i}.pdb"))
        new_pdb = Path(current_dir / (f"{out_prefix}-{i}.pdb"))
        fix_gromacs_pdb(old_pdb, new_pdb, new_chainID=new_chainID)
        old_pdb.unlink()

    if delete_zip_file:
        zip_file.unlink()

    return nframes


def fix_gromacs_pdb(
    pdb_in_fn: Path, pdb_out_fn: Path, *, new_chainID: str = None
) -> None:
    """fix_gromacs_pdb uses Bio to add TER between chains, END at the end,
    and manually make the resSeq numbers continuous.
    Since this changes the PBD numbering it should only be used for structures
    to be used for scoring.
    Specifically, Haddock, that requires chains to have unique resSeq numbers.

    Args:
        pdb_in_fn (Path): Path to input PDB.
        pdb_out_fn (Path): Path to output PDB.
    """

    parsero = PDBParser(QUIET=True)
    pdb = parsero.get_structure("foo", file=pdb_in_fn)

    # Modifying resSeq and chainID
    resSeq = 1
    for resi in pdb.get_residues():
        juan, _, roman = resi._id
        resi._id = (juan, resSeq, roman)
        resSeq += 1
    if new_chainID is not None:
        for chain in pdb.get_chains():
            chain._id = new_chainID

    io = PDBIO()
    io.set_structure(pdb)
    io.save(str(pdb_out_fn))


def join_target_binder(
    pdbs_path: Path,
    nframes: int,
    tar: str = "target",
    bin: str = "binder",
    cpx: str = "complex",
) -> None:
    for i in range(nframes):
        tmp_fn = pdbs_path / f"tmp-{i}.pdb"
        catenate_pdbs(
            tmp_fn,
            FileHandle(pdbs_path / f"{tar}-{i}.pdb"),
            FileHandle(pdbs_path / f"{bin}-{i}.pdb"),
        )
        fix_gromacs_pdb(tmp_fn, pdbs_path / f"{cpx}-{i}.pdb")
        tmp_fn.unlink()


def rm_frames(frames_path: Path, scoring_functions: Iterable, nframes: int) -> None:
    for i in range(nframes):
        Path(frames_path, f"target-{i}.pdb").unlink()
        Path(frames_path, f"binder-{i}.pdb").unlink()
        Path(frames_path, f"complex-{i}.pdb").unlink()

    if "bluues" in scoring_functions:
        bluues_dir = frames_path.parent / "bluues"
        for solv_nrg_file in glob.glob(str(bluues_dir / "*solv_nrg")):
            Path(solv_nrg_file).unlink()
        for pqr_file in glob.glob(str(bluues_dir / "*pqr")):
            Path(pqr_file).unlink()

    if "haddock" in scoring_functions:
        haddock_dir = frames_path.parent / "haddock"
        for pdb_file in glob.glob(str(haddock_dir / "*pdb")):
            Path(pdb_file).unlink()
        for psf_file in glob.glob(str(haddock_dir / "*psf")):
            Path(psf_file).unlink()
