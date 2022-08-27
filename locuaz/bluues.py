from pathlib import Path
from operator import itemgetter
from typing import List
import subprocess as sp
from fileutils import FileHandle, DirHandle
from collections import Sequence
from abstractscoringfunction import AbstractScoringFunction


class bluues(AbstractScoringFunction):
    CPU_TO_MEM_RATIO: int = 4
    bmf_bin_path: FileHandle
    pdb2pqr_bin_path: str = "pdb2pqr30"
    target_chains: tuple
    binder_chains: tuple

    def __init__(
        self, sf_dir, nprocs=2, *, target_chains: Sequence, binder_chains: Sequence
    ):
        self.root_dir = DirHandle(Path(sf_dir, "bluues"), make=False)
        self.nprocs = nprocs
        self.max_concurrent_jobs = self.CPU_TO_MEM_RATIO * self.nprocs

        self.bin_path = FileHandle(self.root_dir / "bluues_new_2")
        self.bmf_bin_path = FileHandle(self.root_dir / "score_bmf_3")

        self.target_chains = tuple(target_chains)
        self.binder_chains = tuple(binder_chains)

    def __call__(self, *, nframes: int, frames_path: Path):

        DirHandle(Path(frames_path, "bluues"), make=True)
        steps = list(range(0, nframes + 1, self.max_concurrent_jobs))
        scores_bluues = []
        scores_bmf = []
        for start, stop in zip(steps[0::1], steps[1::1]):
            step_scores_bluues, step_scores_bmf = self.__submit_batch__(
                start, stop, frames_path
            )
            scores_bluues += step_scores_bluues
            scores_bmf += step_scores_bmf

        # Remaining frames:
        step_scores_bluues, step_scores_bmf = self.__submit_batch__(
            steps[-1], nframes + 1, frames_path
        )
        scores_bluues += step_scores_bluues
        scores_bmf += step_scores_bmf

        return scores_bluues, scores_bmf

    def __submit_batch_pqr__(
        self,
        start: int,
        stop: int,
        frames_path: Path,
    ) -> None:
        processos: List = []
        results_dir = DirHandle(Path(frames_path, "bluues"), make=False)
        for i in range(start, stop):
            pdb_frame = Path(frames_path, "complex-" + str(i) + ".pdb")
            pqr_frame = Path(results_dir, "complex-" + str(i) + ".pqr")

            comando_pdb2pqr = (
                self.pdb2pqr_bin_path
                + " --ff=AMBER "
                + str(pdb_frame)
                + " "
                + str(pqr_frame)
                + " --keep-chain"
            )

            processos.append(
                sp.Popen(
                    comando_pdb2pqr,
                    stdout=sp.PIPE,
                    stderr=sp.PIPE,
                    shell=True,
                    text=True,
                )
            )

        for i, proc in enumerate(processos):
            proc.communicate()

            complex_pqr_frame = Path(results_dir, "complex-" + str(i) + ".pqr")
            target_pqr_frame = Path(results_dir, "target-" + str(i) + ".pqr")
            binder_pqr_frame = Path(results_dir, "binder-" + str(i) + ".pqr")

            f_target = open(target_pqr_frame, "w")
            f_binder = open(binder_pqr_frame, "w")

            with open(complex_pqr_frame, "r") as f_complex:
                for linea in f_complex:
                    if linea[0:4] == "ATOM":
                        if linea[21:22] in self.target_chains:
                            f_target.write(linea)
                        elif linea[21:22] in self.binder_chains:
                            f_binder.write(linea)
            f_target.write("TER\n")
            f_target.write("END")
            f_target.close()

            f_binder.write("TER\n")
            f_binder.write("END")
            f_binder.close()

        return

    def __parse_output__(self, *, score_stdout=None, score_file=None) -> float:
        try:
            with open(str(score_file) + ".solv_nrg", "r") as f:
                for linea in f:
                    if linea[0:26] == "Total              energy:":
                        bluues_raw = float(linea.split()[2])
                        break
        except ValueError as e:
            raise ValueError(f"{self} couldn't parse {score_file}.") from e

        return bluues_raw

    def __get_bluues_bmf_score__(self, results_dir: DirHandle, molecule: str, i: int):
        bluues_file = Path(results_dir, "bluues_" + molecule + "-" + str(i) + ".out")
        bmf_file = Path(results_dir, "bmf_" + molecule + "-" + str(i) + ".out")

        bluues_raw = self.__parse_output__(score_file=bluues_file)
        bluues = bluues_raw * 15 / (abs(bluues_raw) + 15)

        with open(bmf_file, "r") as f:
            lineas = f.readlines()
        bmf, bumps, tors = tuple(map(float, itemgetter(1, 3, 7)(lineas[0].split())))
        gab = 0.17378 * bmf + 0.25789 * bumps + 0.26624 * tors + 0.16446 * bluues

        return bluues_raw, gab

    def __submit_bluues_molecule__(self, results_dir: Path, molecule: str, i: int):
        pqr_frame = Path(results_dir, molecule + "-" + str(i) + ".pqr")
        bluues_output_frame = Path(
            results_dir, "bluues_" + molecule + "-" + str(i) + ".out"
        )
        bmf_output_frame = Path(results_dir, "bmf_" + molecule + "-" + str(i) + ".out")

        comando_bluues = (
            str(self.bin_path) + " " + str(pqr_frame) + " " + str(bluues_output_frame)
        )
        comando_bmf = (
            str(self.bmf_bin_path)
            + " "
            + str(pqr_frame)
            + " "
            + str(bmf_output_frame)
            + " -x"
        )

        proc_bluues = sp.Popen(
            comando_bluues,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            shell=True,
            text=True,
        )
        proc_bmf = sp.Popen(
            comando_bmf,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            shell=True,
            text=True,
        )

        return proc_bluues, proc_bmf

    def __submit_batch__(
        self,
        start: int,
        stop: int,
        frames_path: Path,
    ):
        self.__submit_batch_pqr__(start, stop, frames_path)

        results_dir = DirHandle(Path(frames_path, "bluues"), make=False)
        processos_bluues_cpx: List = []
        processos_bmf_cpx: List = []
        processos_bluues_tar: List = []
        processos_bmf_tar: List = []
        processos_bluues_bin: List = []
        processos_bmf_bin: List = []
        for i in range(start, stop):
            # Complex
            proc_bluues_cpx, proc_bmf_cpx = self.__submit_bluues_molecule__(
                results_dir, "complex", i
            )
            processos_bluues_cpx.append(proc_bluues_cpx)
            processos_bmf_cpx.append(proc_bmf_cpx)

            # Target
            proc_bluues_tar, proc_bmf_tar = self.__submit_bluues_molecule__(
                results_dir, "target", i
            )
            processos_bluues_tar.append(proc_bluues_tar)
            processos_bmf_tar.append(proc_bmf_tar)

            # Binder
            proc_bluues_bin, proc_bmf_bin = self.__submit_bluues_molecule__(
                results_dir, "binder", i
            )
            processos_bluues_bin.append(proc_bluues_bin)
            processos_bmf_bin.append(proc_bmf_bin)

        scores_bluues = []
        scores_bmf = []
        for i, (
            proc_bluues_cpx,
            proc_bmf_cpx,
            proc_bluues_tar,
            proc_bmf_tar,
            proc_bluues_bin,
            proc_bmf_bin,
        ) in enumerate(
            zip(
                processos_bluues_cpx,
                processos_bmf_cpx,
                processos_bluues_tar,
                processos_bmf_tar,
                processos_bluues_bin,
                processos_bmf_bin,
            )
        ):
            proc_bluues_cpx.communicate()
            proc_bmf_cpx.communicate()
            bluues_complex, bmf_complex = self.__get_bluues_bmf_score__(
                results_dir, "complex", i
            )

            proc_bluues_tar.communicate()
            proc_bmf_tar.communicate()
            bluues_target, bmf_target = self.__get_bluues_bmf_score__(
                results_dir, "target", i
            )
            proc_bluues_bin.communicate()
            proc_bluues_bin.communicate()
            bluues_binder, bmf_binder = self.__get_bluues_bmf_score__(
                results_dir, "binder", i
            )

            scores_bluues.append(bluues_complex - bluues_target - bluues_binder)
            scores_bmf.append(bmf_complex - bmf_target - bmf_binder)

        return scores_bluues, scores_bmf
