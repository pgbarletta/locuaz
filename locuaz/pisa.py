from pathlib import Path
from operator import itemgetter
from abc import ABC, ABCMeta
from typing import Any, List
import subprocess as sp
from fileutils import FileHandle, DirHandle
from collections import Sequence
from abstractscoringfunction import AbstractScoringFunction


class pisa(AbstractScoringFunction):

    CPU_TO_MEM_RATIO: int = 1000
    parameters_handle: FileHandle
    target_chains: tuple
    binder_chains: tuple

    def __init__(
        self, sf_dir, nprocs=2, *, target_chains: Sequence, binder_chains: Sequence
    ):
        self.root_dir = DirHandle(Path(sf_dir, "pisa"), make=False)
        self.nprocs = nprocs
        self.max_concurrent_jobs = self.CPU_TO_MEM_RATIO * self.nprocs

        self.parameters_handle = FileHandle(self.root_dir / "pisa.params")
        self.bin_path = FileHandle(self.root_dir / "pisaEnergy_linux")

        self.target_chains = tuple(target_chains)
        self.binder_chains = tuple(binder_chains)

    def __submit_batch__(self, start: int, stop: int, frames_path: Path):
        results_dir = DirHandle(Path(frames_path, "pisa"), make=False)
        processos: List = []

        for i in range(start, stop):
            pdb_frame = Path(frames_path, "complex-" + str(i) + ".pdb")

            comando_pisa = (
                str(self.bin_path)
                + " "
                + str(pdb_frame)
                # fmt: off
                + " " 
                + "".join(self.target_chains) + " " + "".join(self.binder_chains)
                + " "
                # fmt: on
                + str(self.parameters_handle)
            )

            processos.append(
                sp.Popen(
                    comando_pisa,
                    stdout=sp.PIPE,
                    stderr=sp.PIPE,
                    shell=True,
                    text=True,
                )
            )

        scores: List = []
        for proc in processos:
            raw_score, _ = proc.communicate()
            scores.append(float(raw_score))

        return scores

    def __call__(self, *, nframes: int, frames_path: Path):
        print(" -- PISA scoring -- ")
        DirHandle(Path(frames_path, "pisa"), make=True)
        steps = list(range(0, nframes + 1, self.max_concurrent_jobs))
        scores: List = []
        for start, stop in zip(steps[0::1], steps[1::1]):
            scores += self.__submit_batch__(start, stop, frames_path)

        # Remaining frames:
        scores += self.__submit_batch__(steps[-1], nframes + 1, frames_path)

        # clean-up  bach's output file.
        # (Path.cwd() / "output.bss").unlink()

        return scores
