import concurrent.futures as cf
import subprocess as sp
from pathlib import Path
from typing import Tuple, List, Union

from locuaz.abstractscorer import AbstractScorer
from locuaz.complex import GROComplex
from locuaz.fileutils import FileHandle, DirHandle


class Bluues(AbstractScorer):
    bmf_bin_path: FileHandle
    pdb2pqr_bin_path: str = "pdb2pqr30"
    TIMEOUT_PER_FRAME: int = 60

    def __init__(self, sf_dir, *, nthreads=2, mpi_procs=2) -> None:
        super().__init__(sf_dir, nthreads=nthreads, mpi_procs=mpi_procs)
        self.bmf_bin_path = FileHandle(Path(self.root_dir, "bmf"))

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
                    if linea[21:22] == "A":
                        f_target.write(linea)
                    elif linea[21:22] == "B":
                        f_binder.write(linea)
        f_target.write("TER\n")
        f_target.write("END")
        f_target.close()

        f_binder.write("TER\n")
        f_binder.write("END")
        f_binder.close()

        return i

    def __parse_outfile_(self, score_file: Union[Path, FileHandle], original_command: str) -> float:

        with open(Path(score_file), "r") as f:
            for linea in f:
                if linea[0:26] == "Total              energy:":
                    bluues_raw = float(linea.split()[2])
                    return bluues_raw

    def __bluues_molecule__(self, mol: str, i: int) -> float:
        pqr_mol = f"{mol}-{i}.pqr"
        blu_mol_out = f"bluues_{mol}-{i}.out"
        comando_bluues = f"{self.bin_path} {pqr_mol} {blu_mol_out}"
        pbluues = sp.run(
            comando_bluues,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            cwd=self.results_dir,
            shell=True,
            text=True,
        )
        # blu_mol_out_fn = Path(self.results_dir, f"bluues_{mol}-{i}.out")
        blu_mol_out_solv_fn = Path(self.results_dir, f"bluues_{mol}-{i}.out.solv_nrg")
        self.__assert_scorer_outfile__(blu_mol_out_solv_fn, stdout=pbluues.stdout, stderr=pbluues.stderr,
                                                 command=comando_bluues)
        bluues = self.__parse_outfile_(blu_mol_out_solv_fn, comando_bluues)

        return bluues

    def __bluues_worker__(self, i: int) -> Tuple[int, float]:

        bluues_tar = self.__bluues_molecule__("target", i)
        bluues_bin = self.__bluues_molecule__("binder", i)
        bluues_cpx = self.__bluues_molecule__("complex", i)

        bluues = bluues_cpx - bluues_tar - bluues_bin

        return i, bluues

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
        # The first frames will be discarded later.
        scores_bluues: List[float] = [0] * end

        with cf.ProcessPoolExecutor(max_workers=self.nthreads) as exe:
            futuros_pdb2pqr: List[cf.Future] = []
            futuros_bluues: List[cf.Future] = []
            for i in range(start, end):
                futuros_pdb2pqr.append(
                    exe.submit(self.__pdb2pqr_worker__, frames_path, i)
                )
            try:
                timeout = (self.TIMEOUT_PER_FRAME * nframes) // 2
                for futu_pdb2pqr in cf.as_completed(futuros_pdb2pqr, timeout=timeout):
                    if futu_pdb2pqr.exception():
                        print(
                            f"Exception while running pdb2pqr: {futu_pdb2pqr.exception()}",
                            flush=True,
                        )
                        raise futu_pdb2pqr.exception()  # type: ignore

                    j = futu_pdb2pqr.result()
                    futuros_bluues.append(exe.submit(self.__bluues_worker__, j))
            except cf.TimeoutError as e:
                print("pdb2pqr subprocess timed out.", flush=True)
                raise e

            try:
                timeout = self.TIMEOUT_PER_FRAME * nframes
                for futu in cf.as_completed(futuros_bluues, timeout=timeout):
                    if futu.exception():
                        print(
                            f"Exception while running {self.name}: "
                            f"{futu.exception()}",
                            flush=True,
                        )
                        raise futu.exception()  # type: ignore

                    k, bluues = futu.result()
                    scores_bluues[k] = bluues
            except cf.TimeoutError as e:
                print(f"{self.name}/bmf subprocess timed out.", flush=True)
                raise e
        return scores_bluues[start:]
