from asyncio import as_completed
import os
import logging
from pathlib import Path
from socket import timeout
from typing import Tuple, List
import subprocess as sp
from collections.abc import Sequence
import concurrent.futures as cf

from fileutils import FileHandle, DirHandle
from abstractscoringfunction import AbstractScoringFunction


class rosetta(AbstractScoringFunction):
    TIMEOUT_PER_FRAME: int = 30

    def __init__(
        self, sf_dir, nprocs=2, *, target_chains: Sequence, binder_chains: Sequence
    ):
        self.root_dir = DirHandle(Path(sf_dir, "rosetta"), make=False)
        self.bin_path = FileHandle(
            Path(
                self.root_dir,
                "sources/rosetta_source/bin/InterfaceAnalyzer.linuxgccrelease",
            )
        )
        self.executable = (
            str(self.bin_path)
            + " -database "
            + str(Path(self.root_dir, "sources/rosetta_database"))
        )

        self.nprocs = nprocs

        # Set up environment:
        parameters_dir = DirHandle(
            Path(
                self.root_dir,
                "sources/rosetta_source/build/src/release/linux/4.14/64/ppc64le/gcc/8.4",
            ),
            make=False,
        )
        parameters_external_dir = DirHandle(
            Path(
                self.root_dir,
                "sources/rosetta_source/build/external/release/linux/4.14/64/ppc64le/gcc/8.4",
            ),
            make=False,
        )

        os.environ["LD_LIBRARY_PATH"] = (
            os.environ["LD_LIBRARY_PATH"]
            + ":"
            + str(parameters_dir)
            + ":"
            + str(parameters_external_dir)
        )

    def __parse_output__(self, *, score_stdout=None, score_file=None) -> float:
        try:
            with open(score_file, "r") as f:
                lineas = f.readlines()
            score_rosetta = float(lineas[2].split()[3])
        except ValueError as e:
            raise ValueError(f"{self} couldn't parse {score_file}.") from e

        return score_rosetta

    def __rosetta_worker__(self, frames_path: Path, i: int) -> Tuple[int, float]:
        pdb_frame = Path(frames_path, f"complex-{i}.pdb")
        out_rosetta = Path(self.results_dir, f"output_rosetta_{i}.out")

        comando_rosetta = (
            f"{self.executable} --in:file:s {pdb_frame} "
            f"-out:file:score_only {out_rosetta} "
            "-add_regular_scores_to_scorefile true "
            "-overwrite -pack_input true -pack_separated true"
        )

        sp.run(comando_rosetta, stdout=sp.PIPE, stderr=sp.PIPE, shell=True, text=True)
        score_rosetta = self.__parse_output__(score_file=out_rosetta)

        return i, score_rosetta

    def __call__(self, *, nframes: int, frames_path: Path) -> List[float]:

        self.results_dir = DirHandle(Path(frames_path, "rosetta"), make=True)
        scores: List[float] = [0] * (nframes + 1)

        with cf.ProcessPoolExecutor(max_workers=self.nprocs) as exe:
            futuros: List[cf.Future] = []
            for i in range(nframes + 1):
                futuros.append(exe.submit(self.__rosetta_worker__, frames_path, i))

            timeout = self.TIMEOUT_PER_FRAME * nframes
            try:
                for futu in cf.as_completed(futuros, timeout=timeout):
                    if futu.exception():
                        logging.error(
                            f"Exception while running rosetta: {futu.exception()}"
                        )
                        raise futu.exception()  # type: ignore
                    j, score = futu.result()
                    scores[j] = score
            except cf.TimeoutError as e:
                logging.error("rosetta subprocess timed out.")
                raise e
        return scores
