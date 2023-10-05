import os
import concurrent.futures as cf
import subprocess as sp
from pathlib import Path
from typing import Tuple, List

from locuaz.abstractscorer import AbstractScorer
from locuaz.complex import GROComplex
from locuaz.fileutils import FileHandle, DirHandle


class PiePisa(AbstractScorer):
    parameters_handle: FileHandle
    TIMEOUT_PER_FRAME: int = 2

    def __init__(self, sf_dir, *, nthreads=2, mpi_procs=2) -> None:
        self.name = str(self)
        self.root_dir = DirHandle(Path(sf_dir, self.name), make=False)
        self.nthreads = nthreads
        self.mpi_procs = mpi_procs

        self.bin_path_pie = FileHandle(Path(self.root_dir, "pie"))
        self.parameters_handle_pie = FileHandle(Path(self.root_dir, f"pie.params"))
        self.bin_path_pisa = FileHandle(Path(self.root_dir, "pisa"))
        self.parameters_handle_pisa = FileHandle(Path(self.root_dir, f"pisa.params"))

        # Set up environment:
        os.environ["PIEDOCK_HOME"] = str(Path(self.root_dir))

    def __parse_stdout__(self, score_stdout: str, original_command: str) -> float:
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
        comando_pisa = f"{self.bin_path_pisa} {pdb_frame} A B {self.parameters_handle_pisa}"

        p = sp.run(
            comando_pisa,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            cwd=frames_path,
            shell=True,
            text=True,
        )
        pisa_score = self.__parse_stdout__(
            score_stdout=p.stdout, original_command=comando_pisa
        )

        return i, pisa_score

    def __pie_worker__(self, frames_path: Path, i: int) -> Tuple[int, float]:
        pdb_frame = f"complex-{i}.pdb"

        # Using relative path to `pdb_frame` to shorten path to input PDB.
        comando_pie = f"{self.bin_path_pie} {pdb_frame} A B {self.parameters_handle_pie}"

        p = sp.run(
            comando_pie,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            cwd=frames_path,
            shell=True,
            text=True,
        )
        pie_score = self.__parse_stdout__(
            score_stdout=p.stdout, original_command=comando_pie
        )

        return i, pie_score

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
        scores: List[float] = [1] * end

        with cf.ProcessPoolExecutor(max_workers=self.nthreads) as exe:
            futuros: List[cf.Future] = []
            for i in range(start, end):
                futuros.append(exe.submit(self.__pisa_worker__, frames_path, i))
                futuros.append(exe.submit(self.__pie_worker__, frames_path, i))

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
                    scores[j] *= score
            except cf.TimeoutError as e:
                print(f"{self.name} subprocess timed out.", flush=True)
                raise e
            except Exception as e:
                print(f"{self.name} asfdasdfsa.", flush=True)
                raise e
        # Discard the first 0 frames
        return scores[start:]
