from abc import ABCMeta, abstractmethod
from pathlib import Path
from typing import Any, Union, List

from locuaz.complex import GROComplex
from locuaz.fileutils import FileHandle, DirHandle


class AbstractScorer(metaclass=ABCMeta):
    name: str
    root_dir: DirHandle
    results_dir: DirHandle
    bin_path: Union[FileHandle, str]
    nthreads: int
    mpi_procs: int
    TIMEOUT: int = 50

    @abstractmethod
    def __init__(self, sf_dir, *, nthreads=2, mpi_procs=2) -> None:
        self.name = str(self)
        self.root_dir = DirHandle(Path(sf_dir, self.name), make=False)
        self.bin_path = FileHandle(Path(self.root_dir, self.name))
        self.nthreads = nthreads
        self.mpi_procs = mpi_procs

    # Can't annotate output since bluues outputs 2 Lists and the rest 1
    @abstractmethod
    def __call__(self, *, start: int, end: int, frames_path: Path, cpx: GROComplex) -> Any:
        pass

    def __parse_outfile_(self, score_file: Union[Path, FileHandle], original_command: str) -> float:
        pass

    def __parse_outfile_list__(self, score_file: Union[Path, FileHandle], original_command: str) -> List[float]:
        pass

    def __parse_stdout__(self, score_stdout: str, original_command: str) -> float:
        pass

    # Quite hacky, but it works.
    def __str__(self) -> str:
        return str(type(self)).split("'")[1].split(".")[1].lower()

    def __assert_scorer_outfile__(self, score_file: Union[str, Path, FileHandle], *,
                                 stdout: str, stderr: str, command: str) -> Path:
        score_file_path = Path(score_file)
        assert score_file_path.is_file(), f"""{self} error. Can't parse: {score_file}
from:
{command}
with stdout:
{stdout}
and stderr:
{stderr}
"""

        return score_file_path
