from pathlib import Path
from operator import itemgetter
from abc import ABC, ABCMeta
from typing import Any, List
import subprocess as sp
from fileutils import FileHandle, DirHandle
from collections import Sequence
from abstractscoringfunction import AbstractScoringFunction


class evoef2(AbstractScoringFunction):

    CPU_TO_MEM_RATIO: int = 1000
    parameters_handle: FileHandle

    def __init__(
        self, sf_dir, nprocs=2, *, target_chains: Sequence, binder_chains: Sequence
    ):
        self.root_dir = DirHandle(Path(sf_dir, "evoef2"), make=False)
        self.nprocs = nprocs
        self.max_concurrent_jobs = self.CPU_TO_MEM_RATIO * self.nprocs

        self.bin_path = FileHandle(self.root_dir / "EvoEF2")

    def __submit_batch__(self, start: int, stop: int, frames_path: Path):
        results_dir = DirHandle(Path(frames_path, "evoef2"), make=False)
        processos: List = []

        for i in range(start, stop):
            pdb_frame = Path(frames_path, "complex-" + str(i) + ".pdb")

            comando_evoef2 = (
                str(self.bin_path)
                + " --command=ComputeBinding"
                + " --pdb="
                + str(pdb_frame)
            )

            processos.append(
                sp.Popen(
                    comando_evoef2,
                    stdout=sp.PIPE,
                    stderr=sp.PIPE,
                    shell=True,
                    text=True,
                )
            )

        scores: List = []
        for proc in processos:
            raw_score, _ = proc.communicate()
            scores.append(float(raw_score.split()[-4]))

        return scores

    def __call__(self, *, nframes: int, frames_path: Path):
        print(" -- EVOEF2 scoring -- ")
        DirHandle(Path(frames_path, "evoef2"), make=True)
        steps = list(range(0, nframes + 1, self.max_concurrent_jobs))
        scores: List = []
        for start, stop in zip(steps[0::1], steps[1::1]):
            scores += self.__submit_batch__(start, stop, frames_path)

        # Remaining frames:
        scores += self.__submit_batch__(steps[-1], nframes + 1, frames_path)

        # clean-up  bach's output file.
        # (Path.cwd() / "output.bss").unlink()

        return scores
