from pathlib import Path
from typing import Any, List
import subprocess as sp
from fileutils import FileHandle, DirHandle
from abstractscoringfunction import AbstractScoringFunction
from collections import Sequence
import os


class rosetta(AbstractScoringFunction):
    CPU_TO_MEM_RATIO: int = 8

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
        self.max_concurrent_jobs = self.CPU_TO_MEM_RATIO * self.nprocs

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

    def __submit_batch__(
        self,
        start: int,
        stop: int,
        frames_path: Path,
    ):
        processos = []
        for i in range(start, stop):
            pdb_frame = Path(frames_path, "complex-" + str(i) + ".pdb")
            output_rosetta_file = Path(
                self.results_dir, "output_rosetta_" + str(i) + str(".out")
            )

            comando_rosetta = (
                self.executable
                + " --in:file:s "
                + str(pdb_frame)
                + " -out:file:score_only "
                + str(output_rosetta_file)
                + " -add_regular_scores_to_scorefile true"
                + " -overwrite -pack_input true -pack_separated true"
            )

            processos.append(
                sp.Popen(
                    comando_rosetta,
                    stdout=sp.PIPE,
                    stderr=sp.PIPE,
                    shell=True,
                    text=True,
                )
            )

        scores = []
        for i, proc in enumerate(processos):
            sal, err = proc.communicate()
            output_rosetta_file = Path(
                self.results_dir, "output_rosetta_" + str(i) + str(".out")
            )
            score_rosetta = self.__parse_output__(score_file=output_rosetta_file)
            scores.append(score_rosetta)

        return scores

    def __parse_output__(self, *, score_stdout=None, score_file=None) -> float:
        try:
            with open(score_file, "r") as f:
                lineas = f.readlines()
            score_rosetta = float(lineas[2].split()[3])
        except ValueError as e:
            raise ValueError(f"{self} couldn't parse {score_file}.") from e

        return score_rosetta

    def __call__(self, *, nframes: int, frames_path: Path):

        self.results_dir = DirHandle(Path(frames_path, "rosetta"), make=True)
        steps = list(range(0, nframes + 1, self.max_concurrent_jobs))
        scores = []
        for start, stop in zip(steps[0::1], steps[1::1]):
            step_scores = self.__submit_batch__(start, stop, frames_path)
            scores += step_scores

        # Remaining frames:
        step_scores = self.__submit_batch__(steps[-1], nframes + 1, frames_path)
        scores += step_scores

        return scores
