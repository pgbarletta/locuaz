from asyncio import as_completed
import logging
from pathlib import Path
from typing import List, Tuple
import subprocess as sp
from collections.abc import Sequence
import concurrent.futures as cf

from fileutils import FileHandle, DirHandle
from abstractscoringfunction import AbstractScoringFunction


class Evoef2(AbstractScoringFunction):

    CPU_TO_MEM_RATIO: int = 1000
    TIMEOUT_PER_FRAME: int = 1

    def __init__(self, sf_dir, nprocs=2):
        self.root_dir = DirHandle(Path(sf_dir, "evoef2"), make=False)
        self.nprocs = nprocs
        self.max_concurrent_jobs = self.CPU_TO_MEM_RATIO * self.nprocs

        self.bin_path = FileHandle(self.root_dir / "EvoEF2")

    def __evoef2_worker__(self, frames_path: Path, i: int) -> Tuple[int, float]:
        pdb_frame = Path(frames_path, f"complex-{i}.pdb")

        comando_evoef2 = (
            str(self.bin_path)
            + " --command=ComputeBinding"
            + " --pdb="
            + str(pdb_frame.name)
        )

        # Using relative path to `pdb_frame` to shorten path to input PDB.
        p = sp.run(
            comando_evoef2,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            cwd=str(pdb_frame.parent),
            shell=True,
            text=True,
        )

        evoef2_score = self.__parse_output__(score_stdout=p.stdout)

        return i, evoef2_score

    def __parse_output__(self, *, score_stdout=None, score_file=None) -> float:
        try:
            evoef2_score = float(score_stdout.split()[-4])
        except ValueError as e:
            raise ValueError(f"{self} couldn't parse {score_stdout}.") from e

        return evoef2_score

    def __call__(
        self,
        *,
        nframes: int,
        frames_path: Path,
    ) -> List[float]:
        self.results_dir = DirHandle(Path(frames_path, "evoef2"), make=True)
        scores: List[float] = [0] * (nframes)

        with cf.ProcessPoolExecutor(max_workers=self.nprocs) as exe:
            futuros: List[cf.Future] = []

            for i in range(nframes):
                futuros.append(exe.submit(self.__evoef2_worker__, frames_path, i))
            timeout = self.TIMEOUT_PER_FRAME * nframes
            try:
                for futu in cf.as_completed(futuros, timeout=timeout):
                    if futu.exception():
                        print(
                            f"Exception while running evoef2: {futu.exception()}",
                            flush=True,
                        )
                        raise futu.exception()  # type: ignore
                    j, score = futu.result()
                    scores[j] = score
            except cf.TimeoutError as e:
                print("evoef2 subprocess timed out.", flush=True)
                raise e

        return scores
