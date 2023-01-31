from abc import ABCMeta, abstractmethod
from pathlib import Path
from typing import Any, Union, List

from complex import GROComplex
from fileutils import FileHandle, DirHandle


class AbstractScoringFunction(metaclass=ABCMeta):
    name: str
    root_dir: DirHandle
    results_dir: DirHandle
    bin_path: Union[FileHandle, str]
    nprocs: int
    TIMEOUT: int = 50

    @abstractmethod
    def __init__(self, sf_dir, *, nprocs=2) -> None:
        self.name = str(self)
        self.root_dir = DirHandle(Path(sf_dir, self.name), make=False)
        self.bin_path = FileHandle(Path(self.root_dir, self.name))
        self.nprocs = nprocs

    # Can't annotate output since bluues outputs 2 Lists and the rest 1
    @abstractmethod
    def __call__(self, *, nframes: int, frames_path: Path, cpx: GROComplex) -> Any:
        pass

    @abstractmethod
    def __parse_output__(
            self, *, score_stdout: Any = None, score_file: Any = None, original_command=""
    ) -> Union[float, List[float]]:
        pass

    # Quite hacky, but it works.
    def __str__(self) -> str:
        return str(type(self)).split("'")[1].split(".")[1].lower()

    def __assert_scoring_function_outfile__(self, score_file: Union[str, Path, FileHandle], *,
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
