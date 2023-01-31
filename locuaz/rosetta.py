import concurrent.futures as cf
import os
import subprocess as sp
from pathlib import Path
from typing import Tuple, List

from abstractscoringfunction import AbstractScoringFunction
from complex import GROComplex
from fileutils import DirHandle


class Rosetta(AbstractScoringFunction):
    TIMEOUT_PER_FRAME: int = 30

    def __init__(self, sf_dir, *, nprocs=2) -> None:
        super().__init__(sf_dir, nprocs=nprocs)
        self.executable = (
            f'{self.bin_path} -database {Path(self.root_dir, "rosetta_database")}'
        )
        parameters_dir = DirHandle(Path(self.root_dir, "parameters"), make=False)
        parameters_external_dir = DirHandle(
            Path(self.root_dir, "external_parameters"), make=False
        )

        os.environ["LD_LIBRARY_PATH"] = (
                os.environ["LD_LIBRARY_PATH"]
                + ":"
                + str(parameters_dir)
                + ":"
                + str(parameters_external_dir)
        )

    def __parse_output__(
        self, *, score_stdout=None, score_file=None, original_command=""
    ) -> float:

        assert (
            score_file is not None
        ), f"This shouldn't happen. {self} couldn't parse {score_file}\nfrom: \n{original_command}"
        try:
            with open(score_file, "r") as f:
                lineas = f.readlines()
            score_rosetta = float(lineas[2].split()[3])
        except (ValueError, IndexError) as e:
            raise ValueError(
                f"{self} couldn't parse {score_file}\nfrom: \n{original_command}"
            ) from e

        return score_rosetta

    def __rosetta_worker__(self, frames_path: Path, i: int) -> Tuple[int, float]:
        pdb_frame = f"complex-{i}.pdb"
        out_rosetta = f"output_rosetta_{i}.out"
        out_rosetta_fn = Path(self.results_dir, f"output_rosetta_{i}.out")

        comando_rosetta = (
            f"{self.executable} --in:file:s {pdb_frame} "
            f"-out:file:score_only {Path(self.name, out_rosetta)} "
            "-add_regular_scores_to_scorefile true "
            "-overwrite -pack_input true -pack_separated true"
        )

        sp.run(
            comando_rosetta,
            stdout=sp.PIPE,
            cwd=frames_path,
            stderr=sp.PIPE,
            shell=True,
            text=True,
        )
        score_rosetta = self.__parse_output__(
            score_file=out_rosetta_fn, original_command=comando_rosetta
        )

        return i, score_rosetta

    def __call__(
            self,
            *,
            nframes: int,
            frames_path: Path,
            cpx: GROComplex,
    ) -> List[float]:

        self.results_dir = DirHandle(Path(frames_path, self.name), make=True)
        scores: List[float] = [0] * nframes

        with cf.ProcessPoolExecutor(max_workers=self.nprocs) as exe:
            futuros: List[cf.Future] = []
            for i in range(nframes):
                futuros.append(exe.submit(self.__rosetta_worker__, frames_path, i))

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
        return scores
