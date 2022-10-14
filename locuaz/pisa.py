import logging
from pathlib import Path
from typing import Tuple, List, Any
import subprocess as sp
import concurrent.futures as cf

from fileutils import FileHandle, DirHandle
from abstractscoringfunction import AbstractScoringFunction


class Pisa(AbstractScoringFunction):
    name: str = "pisa"
    parameters_handle: FileHandle
    TIMEOUT_PER_FRAME: int = 2

    def __init__(self, sf_dir, nprocs=2):
        self.root_dir = DirHandle(Path(sf_dir, self.name), make=False)
        self.nprocs = nprocs

        self.parameters_handle = FileHandle(Path(self.root_dir, "pisa.params"))
        self.bin_path = FileHandle(Path(self.root_dir, "pisaEnergy_linux"))

    def __parse_output__(
        self, *, score_stdout: Any = None, score_file: Any = None, original_command=""
    ) -> float:
        assert (
            score_stdout is not None
        ), f"This shouldn't happen. {self} couldn't parse {score_stdout}\nfrom: \n{original_command}"

        try:
            pisa_score = float(score_stdout)
        except (ValueError, IndexError) as e:
            raise ValueError(
                f"{self} couldn't parse {score_stdout}\nfrom: \n{original_command}"
            ) from e

        return pisa_score

    def __pisa_worker__(self, frames_path: Path, i: int) -> Tuple[int, float]:
        pdb_frame = f"complex-{i}.pdb"

        # Using relative path to `pdb_frame` to shorten path to input PDB.
        comando_pisa = f"{self.bin_path} {pdb_frame} A B {self.parameters_handle}"

        p = sp.run(
            comando_pisa,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            cwd=frames_path,
            shell=True,
            text=True,
        )
        pisa_score = self.__parse_output__(
            score_stdout=p.stdout, original_command=comando_pisa
        )

        return i, pisa_score

    def __call__(
        self,
        *,
        nframes: int,
        frames_path: Path,
    ) -> List[float]:

        self.results_dir = DirHandle(Path(frames_path, self.name), make=True)
        scores: List[float] = [0] * (nframes)

        with cf.ProcessPoolExecutor(max_workers=self.nprocs) as exe:
            futuros: List[cf.Future] = []
            for i in range(nframes):
                futuros.append(exe.submit(self.__pisa_worker__, frames_path, i))

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
                print("{self.name} subprocess timed out.", flush=True)
                raise e
        return scores
