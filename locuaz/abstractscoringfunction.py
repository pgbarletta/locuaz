from pathlib import Path
from abc import ABCMeta, abstractmethod
from typing import List
from fileutils import FileHandle, DirHandle
from collections.abc import Sequence


class AbstractScoringFunction(metaclass=ABCMeta):
    root_dir: DirHandle
    results_dir: DirHandle
    bin_path: FileHandle
    nprocs: int
    TIMEOUT: int = 50

    @abstractmethod
    def __init__(self, sf_dir, nprocs=2) -> None:
        pass

    # Can't annotate output since bluues outputs 2 Lists and the rest 1
    @abstractmethod
    def __call__(self, *, nframes: int, frames_path: Path):
        pass

    @abstractmethod
    def __parse_output__(self, *, score_stdout=None, score_file=None) -> float:
        pass

    # Quite hacky, but it works.
    def __str__(self) -> str:
        return str(type(self)).split("'")[1].split(".")[1].lower()
