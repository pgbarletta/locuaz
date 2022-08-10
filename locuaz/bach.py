from pathlib import Path
from typing import Any, List
import subprocess as sp
from fileutils import FileHandle, DirHandle
from collections import Sequence
from abstractscoringfunction import AbstractScoringFunction


class bach(AbstractScoringFunction):

    CPU_TO_MEM_RATIO: int = 200
    parameters_handle: FileHandle
    atomic_parameters_handle: FileHandle

    def __init__(
        self, sf_dir, nprocs=2, *, target_chains: Sequence, binder_chains: Sequence
    ):
        self.root_dir = DirHandle(Path(sf_dir, "bach"), make=False)
        self.nprocs = nprocs
        self.max_concurrent_jobs = self.CPU_TO_MEM_RATIO * self.nprocs
        self.parameters_handle = FileHandle(self.root_dir / "BSS.par")
        self.atomic_parameters_handle = FileHandle(
            self.root_dir / "ATOMIC_PARAMETERS_BSS"
        )
        self.bin_path = FileHandle(self.root_dir / "bach")

    def __submit_batch__(self, start: int, stop: int, frames_path: Path):
        results_dir = DirHandle(Path(frames_path, "bach"), make=False)
        processos: List = []

        for i in range(start, stop):
            pdb_frame = Path(frames_path, "complex-" + str(i) + ".pdb")
            pdb_frame_txt = Path(str(results_dir), "pdb_to_bach_" + str(i) + ".txt")
            with open(pdb_frame_txt, "w") as f:
                f.write(str(pdb_frame))

            comando_bach = (
                str(self.bin_path)
                + " -COMPUTE_ENE -STRICT_INTERFACE -PDBLIST "
                + str(pdb_frame_txt)
                + " -FILE_PAR "
                + str(self.parameters_handle)
                + " -FILE_PAR_AT "
                + str(self.atomic_parameters_handle)
            )

            processos.append(
                sp.Popen(
                    comando_bach,
                    stdout=sp.PIPE,
                    stderr=sp.PIPE,
                    shell=True,
                    text=True,
                )
            )

        scores_bach: List = []
        for proc in processos:
            raw_score, _ = proc.communicate()
            scores_bach.append(float(raw_score.split()[-1]))

        return scores_bach

    def __call__(self, *, nframes: int, frames_path: Path):
        print(" -- BACH scoring -- ")
        DirHandle(Path(frames_path, "bach"), make=True)
        steps = list(range(0, nframes + 1, self.max_concurrent_jobs))
        scores: List = []
        for start, stop in zip(steps[0::1], steps[1::1]):
            scores += self.__submit_batch__(start, stop, frames_path)

        # Remaining frames:
        scores += self.__submit_batch__(steps[-1], nframes + 1, frames_path)

        # clean-up  bach's output file.
        (Path.cwd() / "output.bss").unlink()

        return scores

    def __str__(self) -> str:
        return "bach"
