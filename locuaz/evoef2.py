from pathlib import Path
from typing import List
import subprocess as sp
from fileutils import FileHandle, DirHandle
from collections import Sequence
from abstractscoringfunction import AbstractScoringFunction


class evoef2(AbstractScoringFunction):

    CPU_TO_MEM_RATIO: int = 1000

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
            evoef2_score = self.__parse_output__(score_stdout=raw_score)
            scores.append(evoef2_score)

        return scores

    def __parse_output__(self, *, score_stdout=None, score_file=None) -> float:
        try:
            evoef2_score = float(score_stdout.split()[-4])
        except ValueError as e:
            raise ValueError(f"{self} couldn't parse {score_stdout}.") from e

        return evoef2_score

    def __call__(self, *, nframes: int, frames_path: Path):
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
