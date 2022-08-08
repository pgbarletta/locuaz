from pathlib import Path
from collections import abc
from molecules import GROComplex, PDBStructure
from biobb_md.gromacs.gmxselect import Gmxselect
from fileutils import FileHandle, update_header, catenate


def update_all_ndxs(cpx: GROComplex):
    try:
        cpx.target_ndx = write_chain_selection_ndx(
            "target",
            cpx.top.target_chains.keys(),
            cpx.pdb,
            cpx.dir_handle.dir_path,
        )
        cpx.binder_ndx = write_chain_selection_ndx(
            "binder",
            cpx.top.binder_chains.keys(),
            cpx.pdb,
            cpx.dir_handle.dir_path,
        )

        cpx.complex_ndx = catenate(
            cpx.dir_handle.dir_path / "complex.ndx", cpx.target_ndx, cpx.binder_ndx
        )
    except Exception as e:
        raise e


def write_chain_selection_ndx(
    selname: str, chains: abc.Sequence, pdb: PDBStructure, out_dir: Path
) -> FileHandle:

    selection_target = " or ".join([f"chain {chainID}" for chainID in chains])
    prop = {"selection": selection_target}
    ndx_fn = Path(out_dir) / f"{selname}.ndx"

    comando = Gmxselect(
        input_structure_path=str(pdb.file.path),
        output_ndx_path=str(ndx_fn),
        properties=prop,
    )
    comando.can_write_console_log = False
    Error_code = comando.launch()
    if Error_code != 0:
        print(
            f"Could not write .ndx files. Error_code: {Error_code}. Check log file.",
            flush=True,
        )
        raise RuntimeError
    ndx = FileHandle(ndx_fn)
    update_header(ndx, f"[ {selname} ]\n")

    return ndx


def write_non_overlapping_ndx(
    pdb_path: Path, ndx_fn_in: FileHandle, resid: int, dist_threshold: float = 0.3
) -> FileHandle:

    # First, get the number of waters that overlap
    options_1 = {
        "selection": f" (same residue as resname SOL and "
        f"within {dist_threshold} of (group 0 and resid {resid}))"
    }
    ndx_fn_wat_out = str(pdb_path.parent / "overlapping.ndx")
    select_overlapping_waters = Gmxselect(
        input_structure_path=str(pdb_path),
        input_ndx_path=str(ndx_fn_in),
        output_ndx_path=str(ndx_fn_wat_out),
        properties=options_1,
    )
    select_overlapping_waters.can_write_console_log = False
    Error_code = select_overlapping_waters.launch()
    if Error_code != 0:
        print(
            f"Could not write .ndx files. Error_code: {Error_code}. Check log file.",
            flush=True,
        )
        raise RuntimeError

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
    options_1 = {
        "selection": f"not (same residue as resname SOL and "
        f"within {dist_threshold} of (group 0 and resid {resid}))"
    }
    ndx_fn_out = str(pdb_path.parent / "non_overlapping.ndx")
    select_overlapping_waters = Gmxselect(
        input_structure_path=str(pdb_path),
        input_ndx_path=str(ndx_fn_in),
        output_ndx_path=str(ndx_fn_out),
        properties=options_1,
    )
    select_overlapping_waters.can_write_console_log = False
    Error_code = select_overlapping_waters.launch()
    if Error_code != 0:
        print(
            f"Could not write .ndx files. Error_code: {Error_code}. Check log file.",
            flush=True,
        )
        raise RuntimeError
    ndx = FileHandle(ndx_fn_out)
    update_header(ndx, f"[ non_overlapping ]\n")

    return ndx, wat_atoms // 3
