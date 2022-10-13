import logging
from pathlib import Path
from typing import Tuple, List
import subprocess as sp
from collections.abc import Sequence
import concurrent.futures as cf

from fileutils import FileHandle, DirHandle
from abstractscoringfunction import AbstractScoringFunction


class Bach(AbstractScoringFunction):

    name: str = "bach"
    parameters_handle: FileHandle
    atomic_parameters_handle: FileHandle
    TIMEOUT_PER_FRAME: int = 2

    def __init__(self, sf_dir, nprocs=2):
        self.root_dir = DirHandle(Path(sf_dir, self.name), make=False)
        self.nprocs = nprocs
        self.parameters_handle = FileHandle(self.root_dir / "BSS.par")
        self.atomic_parameters_handle = FileHandle(
            self.root_dir / "ATOMIC_PARAMETERS_BSS"
        )
        self.bin_path = FileHandle(self.root_dir / self.name)

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
        bach_score = self.__parse_output__(
            score_stdout=p.stdout, original_command=comando_bach
        )

        return i, bach_score

    def __parse_output__(
        self, *, score_stdout=None, score_file=None, original_command=""
    ) -> float:
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
        nframes: int,
        frames_path: Path,
    ) -> List[float]:

        self.results_dir = DirHandle(Path(frames_path, self.name), make=True)
        scores: List[float] = [0] * (nframes)

        with cf.ProcessPoolExecutor(max_workers=self.nprocs) as exe:
            futuros: List[cf.Future] = []
            for i in range(nframes):
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
                print("bach subprocess timed out.", flush=True)
                raise e
        try:
            (Path.cwd() / "output.bss").unlink()
        except FileNotFoundError:
            pass

        return scores
