import concurrent.futures as cf
import subprocess as sp
from pathlib import Path
from typing import List, Tuple, Union
import re

from locuaz.abstractscorer import AbstractScorer
from locuaz.complex import GROComplex
from locuaz.fileutils import DirHandle, FileHandle


class AutodockVina(AbstractScorer):
    CPU_TO_MEM_RATIO: int = 1000
    TIMEOUT_PER_FRAME: int = 1

    def __init__(self, sf_dir, *, nthreads=2, mpi_procs=2) -> None:
        super().__init__(sf_dir, nthreads=nthreads, mpi_procs=mpi_procs)
        self.max_concurrent_jobs = self.CPU_TO_MEM_RATIO * self.nthreads
        self.openbabel_bin = "obabel"
        self.rgx = re.compile("Estimated Free Energy of Binding\s*:\s([+-]?([0-9]*[.])?[0-9]+)")

    def __autodockvina_worker__(self, i: int) -> Tuple[int, float]:
        # Get PDBQT file for target
        target_pdb = f"../target-{i}.pdb"
        target_pdbqt = f"target-{i}.pdbqt"

        comando_ob_target = \
            f"{self.openbabel_bin} -ipdb {target_pdb} -O {target_pdbqt} --partialcharge gasteiger  ---errorlevel 0"

        p = sp.run(comando_ob_target, stdout=sp.PIPE, stderr=sp.PIPE, cwd=self.results_dir, shell=True,text=True)

        self.__assert_scorer_outfile__(Path(self.results_dir, target_pdbqt), stdout=p.stdout, stderr=p.stderr,
                                                 command=comando_ob_target)
        # Get PDBQT file for binder
        binder_pdb = f"../binder-{i}.pdb"
        binder_pdbqt = f"binder-{i}.pdbqt"

        comando_ob_binder = \
            f"{self.openbabel_bin} -ipdb {binder_pdb} -O {binder_pdbqt} --partialcharge gasteiger -xr ---errorlevel 0"

        p = sp.run(comando_ob_binder, stdout=sp.PIPE, stderr=sp.PIPE, cwd=self.results_dir, shell=True, text=True)

        self.__assert_scorer_outfile__(Path(self.results_dir, binder_pdbqt), stdout=p.stdout, stderr=p.stderr,
                                                 command=comando_ob_binder)

        comando_vina = f"{self.bin_path} --score_only --ligand {target_pdbqt} --receptor {binder_pdbqt} --autobox"

        p = sp.run(comando_vina, stdout=sp.PIPE, stderr=sp.PIPE, cwd=self.results_dir, shell=True, text=True)

        autodockvina_score = self.__parse_stdout__(
            score_stdout=p.stdout, original_command=comando_vina
        )

        return i, autodockvina_score

    def __parse_stdout__(self, score_stdout: str, original_command: str) -> float:

        try:
            autodockvina_score = float(self.rgx.search(score_stdout).group(1))
        except (ValueError, IndexError, Exception) as e:
            raise ValueError(
                f"{self} couldn't parse {score_stdout}\nfrom: \n{original_command}"
            ) from e

        return autodockvina_score

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
                futuros.append(exe.submit(self.__autodockvina_worker__, i))
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

        # Discard the first 0 frames
        return scores[start:]

    def __parse_outfile_(self, score_file: Union[Path, FileHandle], original_command: str) -> float:
        pass
