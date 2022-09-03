from pathlib import Path
from operator import itemgetter
from typing import Tuple, List
import subprocess as sp
import logging
from collections.abc import Sequence
import concurrent.futures as cf

from fileutils import FileHandle, DirHandle
from abstractscoringfunction import AbstractScoringFunction


class bluues(AbstractScoringFunction):
    bmf_bin_path: FileHandle
    pdb2pqr_bin_path: str = "pdb2pqr30"
    target_chains: tuple
    binder_chains: tuple
    TIMEOUT_PER_FRAME: int = 60

    def __init__(
        self, sf_dir, nprocs=2, *, target_chains: Sequence, binder_chains: Sequence
    ):
        self.root_dir = DirHandle(Path(sf_dir, "bluues"), make=False)
        self.nprocs = nprocs

        self.bin_path = FileHandle(self.root_dir / "bluues_new_2")
        self.bmf_bin_path = FileHandle(self.root_dir / "score_bmf_3")

        self.target_chains = tuple(target_chains)
        self.binder_chains = tuple(binder_chains)

    def __pdb2pqr_worker__(self, frames_path: Path, i: int) -> int:

        pdb_frame = Path(frames_path, "complex-" + str(i) + ".pdb")
        pqr_frame = Path(self.results_dir, "complex-" + str(i) + ".pqr")

        comando_pdb2pqr = (
            self.pdb2pqr_bin_path
            + " --ff=AMBER "
            + str(pdb_frame)
            + " "
            + str(pqr_frame)
            + " --keep-chain"
        )
        sp.run(comando_pdb2pqr, stdout=sp.PIPE, stderr=sp.PIPE, shell=True, text=True)
        complex_pqr_frame = Path(self.results_dir, "complex-" + str(i) + ".pqr")
        target_pqr_frame = Path(self.results_dir, "target-" + str(i) + ".pqr")
        binder_pqr_frame = Path(self.results_dir, "binder-" + str(i) + ".pqr")

        # Split the complex .pqr file into target and binder .pqr files.
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

        return i

    def __parse_output__(
        self, *, score_stdout=None, score_file=None
    ) -> Tuple[float, float]:
        try:
            with open(str(score_file[0]) + ".solv_nrg", "r") as f:
                for linea in f:
                    if linea[0:26] == "Total              energy:":
                        bluues_raw = float(linea.split()[2])
                        break
        except ValueError as e:
            raise ValueError(f"{self} couldn't parse {score_file[0]}.") from e
        bluues = bluues_raw * 15 / (abs(bluues_raw) + 15)

        try:
            with open(score_file[1], "r") as f:
                lineas = f.readlines()
            bmf, bumps, tors = tuple(map(float, itemgetter(1, 3, 7)(lineas[0].split())))
            gab = 0.17378 * bmf + 0.25789 * bumps + 0.26624 * tors + 0.16446 * bluues
        except ValueError as e:
            raise ValueError(f"{self} couldn't parse {score_file[1]}.") from e

        return bluues_raw, gab

    def __bluues_bmf_molecule__(self, mol: str, i: int) -> Tuple[float, float]:
        pqr_mol = Path(self.results_dir, f"{mol}-{i}.pqr")
        blu_mol_out = Path(self.results_dir, f"bluues_{mol}-{i}.out")

        # BLUUES
        comando_bluues = f"{self.bin_path} {pqr_mol} {blu_mol_out}"
        sp.run(
            comando_bluues,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            shell=True,
            text=True,
        )

        # BMF
        bmf_mol_out = Path(self.results_dir, f"bmf_{mol}-{i}.out")
        comando_bmf = f"{self.bmf_bin_path} {pqr_mol} {bmf_mol_out} -x"
        sp.run(
            comando_bmf,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            shell=True,
            text=True,
        )
        bluues, bmf = self.__parse_output__(score_file=(blu_mol_out, bmf_mol_out))

        return bluues, bmf

    def __bluues_bmf_worker__(self, i: int) -> Tuple[int, float, float]:

        bluues_tar, bmf_tar = self.__bluues_bmf_molecule__("target", i)
        bluues_bin, bmf_bin = self.__bluues_bmf_molecule__("binder", i)
        bluues_cpx, bmf_cpx = self.__bluues_bmf_molecule__("complex", i)

        bluues = bluues_cpx - bluues_tar - bluues_bin
        bmf = bmf_cpx - bmf_tar - bmf_bin

        return i, bluues, bmf

    def __call__(
        self, *, nframes: int, frames_path: Path
    ) -> Tuple[List[float], List[float]]:

        self.results_dir = DirHandle(Path(frames_path, "bluues"), make=True)
        scores_bluues: List[float] = [0] * (nframes + 1)
        scores_bmf: List[float] = [0] * (nframes + 1)
        # TODO: check bluues doesn't do anything weird.
        with cf.ProcessPoolExecutor(max_workers=self.nprocs) as exe:
            futuros_pdb2pqr: List[cf.Future] = []
            futuros_bluues_bmf: List[cf.Future] = []
            for i in range(nframes + 1):
                futuros_pdb2pqr.append(
                    exe.submit(self.__pdb2pqr_worker__, frames_path, i)
                )
            try:
                timeout = (self.TIMEOUT_PER_FRAME * nframes) // 2
                for futu_pdb2pqr in cf.as_completed(futuros_pdb2pqr, timeout=timeout):
                    if futu_pdb2pqr.exception():
                        logging.error(
                            f"Exception while running pdb2pqr: {futu_pdb2pqr.exception()}"
                        )
                        raise futu_min.exception()  # type: ignore

                    j = futu_pdb2pqr.result()
                    futuros_bluues_bmf.append(exe.submit(self.__bluues_bmf_worker__, j))
            except cf.TimeoutError as e:
                logging.error("pdb2pqr subprocess timed out.")
                raise e

            try:
                timeout = self.TIMEOUT_PER_FRAME * nframes
                for futu in cf.as_completed(futuros_bluues_bmf, timeout=timeout):
                    if futu.exception():
                        logging.error(
                            f"Exception while running pdb2pqr: " f"{futu.exception()}"
                        )
                        raise futu.exception()  # type: ignore

                    k, bluues, bmf = futu.result()
                    scores_bluues[k] = bluues
                    scores_bmf[k] = bmf
            except cf.TimeoutError as e:
                logging.error("bluues/bmf subprocess timed out.")
                raise e
        return scores_bluues, scores_bmf
