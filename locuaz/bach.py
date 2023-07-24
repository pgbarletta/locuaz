import concurrent.futures as cf
import subprocess as sp
from pathlib import Path
from typing import Tuple, List

from locuaz.abstractscorer import AbstractScorer
from locuaz.complex import GROComplex
from locuaz.fileutils import FileHandle, DirHandle


class Bach(AbstractScorer):
    parameters_handle: FileHandle
    atomic_parameters_handle: FileHandle
    TIMEOUT_PER_FRAME: int = 2

    def __init__(self, sf_dir, *, nthreads=2, mpi_procs=2) -> None:
        super().__init__(sf_dir, nthreads=nthreads, mpi_procs=mpi_procs)
        self.parameters_handle = FileHandle(Path(self.root_dir, "BSS.par"))
        self.atomic_parameters_handle = FileHandle(
            Path(self.root_dir, "ATOMIC_PARAMETERS_BSS")
        )

    def __bach_worker__(self, frames_path: Path, i: int) -> Tuple[int, float]:
        # Use relative paths to shorten the input.
        pdb_frame = f"complex-{i}.pdb"
        pdb_frame_txt = f"pdb_to_bach_{i}.txt"
        with open(Path(self.results_dir, pdb_frame_txt), "w") as f:
            f.write(str(pdb_frame))

        comando_bach = (
            str(self.bin_path)
            + " -COMPUTE_ENE -STRICT_INTERFACE -PDBLIST "
            + str(Path(self.name, pdb_frame_txt))
            + " -FILE_PAR "
            + str(self.parameters_handle)
            + " -FILE_PAR_AT "
            + str(self.atomic_parameters_handle)
        )
        p = sp.run(
            comando_bach,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            cwd=frames_path,
            shell=True,
            text=True,
        )
        bach_score = self.__parse_stdout__(
            score_stdout=p.stdout, original_command=comando_bach
        )

        return i, bach_score

    def __parse_stdout__(self, score_stdout: str, original_command: str) -> float:
        assert (
            score_stdout is not None
        ), f"This shouldn't happen. {self} couldn't parse {score_stdout}\nfrom: \n{original_command}"

        try:
            bach_score = float(score_stdout.split()[-1])
        except (ValueError, IndexError) as e:
            raise ValueError(
                f"{self} couldn't parse {score_stdout}\nfrom: \n{original_command}"
            ) from e

        return bach_score

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
                futuros.append(exe.submit(self.__bach_worker__, frames_path, i))

            timeout = self.TIMEOUT_PER_FRAME * nframes
            try:
                for futu in cf.as_completed(futuros, timeout=timeout):
                    if futu.exception():
                        print(
                            f"Exception while running bach: {futu.exception()}",
                            flush=True,
                        )
                        raise futu.exception()  # type: ignore
                    j, score = futu.result()
                    scores[j] = score
            except cf.TimeoutError as e:
                print(f"bach subprocess timed out.", flush=True)
                raise e
        try:
            (Path.cwd() / "output.bss").unlink()
        except FileNotFoundError:
            pass

        # Discard the first 0 frames
        return scores[start:]
