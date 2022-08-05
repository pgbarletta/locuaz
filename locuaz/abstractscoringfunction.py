from pathlib import Path
from abc import ABCMeta
from typing import List
from fileutils import FileHandle, DirHandle
from collections import Sequence


class AbstractScoringFunction(metaclass=ABCMeta):
    root_dir: DirHandle
    results_dir: DirHandle
    bin_path: FileHandle

    nprocs: int
    # Some scoring fuctions (eg: bluues) fill too much RAM, so the number of
    # concurrent jobs, have to be limited accordingly:
    max_concurrent_jobs: int
    CPU_TO_MEM_RATIO: int = 5

    def __init__(
        self, sf_dir, nprocs=2, *, target_chains: Sequence, binder_chains: Sequence
    ):
        pass

    def __call__(self, *, nframes: int, frames_path: Path) -> List:
        pass
