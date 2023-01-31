import csv
import shutil as sh
import subprocess as sp
import zipfile
from pathlib import Path
from typing import List, Any

from abstractscoringfunction import AbstractScoringFunction
from complex import GROComplex
from fileutils import DirHandle, FileHandle
from molecules import ZipTopology


class Gmx_mmpbsa(AbstractScoringFunction):
    bin_name = "gmx_MMPBSA"
    TIMEOUT_PER_FRAME: int = 20

    def __init__(self, sf_dir, *, nprocs=2) -> None:
        super().__init__(sf_dir, nprocs=nprocs)
        # `gmx_mmpbsa` isn't actually the binary, but the config file
        self.in_path = self.bin_path
        if nprocs > 1:
            self.bin_name = f"mpirun -np {min(48, nprocs)} gmx_MMPBSA MPI"

    def __parse_output__(
        self, *, score_stdout: Any = None, score_file: Any = None, original_command=""
    ) -> List[float]:

        with open(Path(score_file), 'r') as csv_file:
            text = csv.reader(csv_file)
            for line in text:
                if line == ["Delta Energy Terms"]:
                    next(text)
                    break
                else:
                    continue
            mmpbsa_score: List[float] = [float(line[-1]) for line in text if len(line) > 0]
        return mmpbsa_score

    def __call__(
        self,
        *,
        nframes: int,
        frames_path: Path,
        cpx: GROComplex,
    ) -> List[float]:

        results_dir = self.__initialize_scoring_dir__(frames_path, cpx)
        score_gmxmmpbsa = Path(results_dir, "score_gmxmmpbsa.csv")
        comando_gmx_MMPBSA = f"{self.bin_name} -O -i {self.in_path} -cp {self.top} -cs {cpx.tpr} "\
                             f"-ci {cpx.ndx} -cg 0 1 -ct {self.trj} -eo {score_gmxmmpbsa} -nogui"

        try:
            p = sp.run(
                comando_gmx_MMPBSA,
                stdout=sp.PIPE,
                stderr=sp.PIPE,
                cwd=Path(results_dir),
                shell=True,
                text=True,
                timeout=self.TIMEOUT_PER_FRAME * nframes,
            )
        except RuntimeError as e:
            print(f"{self.name} subprocess timed out.", flush=True)
            raise e

        self.__assert_scoring_function_outfile__(score_gmxmmpbsa, stdout=p.stdout, stderr=p.stderr,
                                                 command=comando_gmx_MMPBSA)
        mmpbsa_score = self.__parse_output__(score_file=score_gmxmmpbsa, original_command=comando_gmx_MMPBSA)

        return mmpbsa_score

    def __initialize_scoring_dir__(self, frames_path: Path, cpx: GROComplex,) -> DirHandle:
        results_dir = DirHandle(Path(frames_path, self.name), make=True)

        if isinstance(cpx.top, ZipTopology):
            zipped_top = zipfile.ZipFile(Path(cpx.top))
            try:
                top_name = next(filter(lambda x: x[-4:] == '.top', zipped_top.namelist()))
            except StopIteration:
                raise RuntimeError(f"Could not find a .top file in {cpx.top}")
            with zipped_top as top:
                top.extractall(results_dir)
            self.top = FileHandle(Path(results_dir, top_name))
        else:
            print("No ZipTopology, are you sure?")
            sh.copy(cpx.top, Path(results_dir))
            self.top = cpx.top.file
        try:
            ext = Path(cpx.tra).suffix
            self.trj = FileHandle(Path(frames_path, f"fix_{cpx.name}{ext}"))
        except FileNotFoundError as e:
            raise RuntimeError(f"gmx_mmpbsa: cannot find fix_{cpx.name} in {frames_path}") from e

        return results_dir
