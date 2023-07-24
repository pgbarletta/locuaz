import concurrent.futures as cf
import subprocess as sp
from pathlib import Path
from typing import List, Tuple

from locuaz.abstractscorer import AbstractScorer
from locuaz.complex import GROComplex
from locuaz.fileutils import DirHandle


class Evoef2(AbstractScorer):
    CPU_TO_MEM_RATIO: int = 1000
    TIMEOUT_PER_FRAME: int = 1

    def __init__(self, sf_dir, *, nthreads=2, mpi_procs=2) -> None:
        super().__init__(sf_dir, nthreads=nthreads, mpi_procs=mpi_procs)
        self.max_concurrent_jobs = self.CPU_TO_MEM_RATIO * self.nthreads

    def __evoef2_worker__(self, frames_path: Path, i: int) -> Tuple[int, float]:
        # Using relative path to `pdb_frame` to shorten path to input PDB.
        pdb_frame = f"complex-{i}.pdb"

        comando_evoef2 = f"{self.bin_path} --command=ComputeBinding --pdb={pdb_frame}"

        p = sp.run(
            comando_evoef2,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            cwd=frames_path,
            shell=True,
            text=True,
        )

        evoef2_score = self.__parse_stdout__(
            score_stdout=p.stdout, original_command=comando_evoef2
        )

        return i, evoef2_score

    def __parse_stdout__(self, score_stdout: str, original_command: str) -> float:
        assert (
            score_stdout is not None
        ), f"This shouldn't happen. {self} couldn't parse {score_stdout}\nfrom: \n{original_command}"

        try:
            evoef2_score = float(score_stdout.split()[-4])
        except (ValueError, IndexError) as e:
            raise ValueError(
                f"{self} couldn't parse {score_stdout}\nfrom: \n{original_command}"
            ) from e

        return evoef2_score

    def __call__(
            self,
            *,
            start: int,
            end: int,
            frames_path: Path,
            cpx: GROComplex,
    ) -> List[float]:
        self.results_dir = DirHandle(Path(frames_path, self.name), make=True)
        nframes = end - start
        # The first unused frames will be discarded later
        scores: List[float] = [0] * end

        with cf.ProcessPoolExecutor(max_workers=self.nthreads) as exe:
            futuros: List[cf.Future] = []

            for i in range(start, end):
                futuros.append(exe.submit(self.__evoef2_worker__, frames_path, i))
            timeout = self.TIMEOUT_PER_FRAME * nframes
            try:
                for futu in cf.as_completed(futuros, timeout=timeout):
                    if futu.exception():
                        print(
                            f"Exception while running {self.name}: {futu.exception()}",
                            flush=True,
                        )
                        raise futu.exception()  # type: ignore
                    j, score = futu.result()
                    scores[j] = score
            except cf.TimeoutError as e:
                print(f"{self.name} subprocess timed out.", flush=True)
                raise e

        # Discard the first 0 frames
        return scores[start:]
