from typing import Dict, Optional, Generic, TypeVar, Tuple, List, Iterator
from collections.abc import Sequence, Iterable
from pathlib import Path
import bisect
import shutil as sh
from warnings import warn

from Bio.SeqUtils import seq1

# This will be used to map non-conventional AAs to conventional ones, so the
# scorers.
AA_MAP: Dict[str, str] = {
    "ALA": "ALA",
    "ARG": "ARG",
    "ASN": "ASN",
    "ASP": "ASP",
    "CYS": "CYS",
    "CYX": "CYS",
    "CY2": "CYS",
    "GLN": "GLN",
    "GLU": "GLU",
    "GLY": "GLY",
    "HIS": "HIS",
    "HIE": "HIS",
    "HE1": "HIS",
    "HID": "HIS",
    "ILE": "ILE",
    "LEU": "LEU",
    "LYS": "LYS",
    "MET": "MET",
    "PHE": "PHE",
    "PRO": "PRO",
    "SER": "SER",
    "THR": "THR",
    "TRP": "TRP",
    "TYR": "TYR",
    "VAL": "VAL",
}

T = TypeVar("T")


class PriorityDeque(Generic[T], Sequence, Iterable):
    """
    PriorityDeque is a fixed-size Priority queue, unless the input maxsize is negative
    in which case it has not size limit.
    I tried to implement it using heappush, heappop and heappushpop, but it uses
    a weird ordering, and it forces the element type to be comparable as well,
    to decide between pairs that have the same priorities. In this class, that's
    left to bisect, since we don't care about sorting stability.
    """

    def __init__(self, maxsize: int):
        self.maxsize: int = maxsize
        self.priorities: List[int] = []
        self.elements: List[T] = []

    def put(self, pair: Tuple[int, T]) -> None:
        try:
            idx = bisect.bisect(self.priorities, int(pair[0]))
        except ValueError:
            raise (
                f"{pair[0]} is a {type(pair[0])}. Input must be a tuple with"
                "an int as first element."
            )
        self.priorities.insert(idx, pair[0])
        self.elements.insert(idx, pair[1])
        if len(self.elements) > self.maxsize:
            self.priorities.pop()
            self.elements.pop()

    def get(self) -> Tuple[int, T]:
        return (self.priorities.pop(0), self.elements.pop(0))

    def __getitem__(self, item) -> Tuple[int, T]:
        return (self.priorities[item], self.elements[item])

    def __iter__(self) -> Iterator:
        return zip(self.priorities, self.elements)

    def __len__(self):
        return len(self.elements)

    def __str__(self) -> str:
        return str(list(zip(self.priorities, self.elements)))


class GromacsError(Exception):
    pass


class MutationError(Exception):
    pass


class GenerationError(Exception):
    pass


class UserInputError(Exception):
    pass


def launch_biobb(
    biobb_obj,
    *,
    can_write_console_log: bool = False,
    backup_dict: Optional[Path] = None,
) -> None:
    """

    Parameters
    ----------
    biobb_obj
    can_write_console_log
    backup_dict
    """
    biobb_obj.can_write_console_log = can_write_console_log
    err = biobb_obj.launch()
    if err == 0 or err is None:
        try:
            Path(biobb_obj.out_log.name).unlink()
            Path(biobb_obj.err_log.name).unlink()
        except (FileNotFoundError, Exception):
            # log files may get the same name when running parallel biobb processes
            # and one run may try to delete a log file that was already deleted.
            pass
    else:
        log_out = Path(biobb_obj.out_log.name)
        log_err = Path(biobb_obj.err_log.name)
        if backup_dict:
            try:
                new_log_out = backup_dict / f"{log_out.name}_out.txt"
                sh.move(log_out, new_log_out)
                log_out = new_log_out

                new_log_err = backup_dict / f"{log_err.name}_err.txt"
                sh.move(log_err, new_log_err)
                log_err = new_log_err
            except FileNotFoundError:
                # There may be concurrent biobb processes so this file may have been copied already. In this case,
                # the user will probably get a log file in the wrong `backup_dict`, but I can't do much about this
                pass

        raise GromacsError(f"{biobb_obj} failed. Check {log_out} and {log_err} .")


def ext(name: str, suffix: str) -> str:
    """ext utility function, so I don't have to worry about the extension.

    Args:
        name (str): filename with or without the extension.
        suffix (str): desired extension.

    Returns:
        str: filename with the extension.
    """
    return f"{name.split('.')[0]}.{suffix}"


def my_seq1(resn_3: str) -> str:
    resn_1 = seq1(resn_3)
    if resn_1 == "X":
        try:
            resn_1 = seq1(AA_MAP[resn_3])
            warn(
                f"Converted non-standard residue {resn_3} from 3-letter code to 1-letter: {resn_1}."
            )
        except KeyError:
            warn(
                f"Could not convert non-standard residue {resn_3} from 3-letter code to 1-letter. Setting it to 'X'."
            )
    return resn_1
