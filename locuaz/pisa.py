import logging
from pathlib import Path
from typing import Tuple, List
import subprocess as sp
from collections.abc import Sequence
import concurrent.futures as cf

from fileutils import FileHandle, DirHandle
from abstractscoringfunction import AbstractScoringFunction


class pisa(AbstractScoringFunction):

    parameters_handle: FileHandle
    target_chains: tuple
    binder_chains: tuple
    TIMEOUT_PER_FRAME: int = 2

    def __init__(
        self, sf_dir, nprocs=2, *, target_chains: Sequence, binder_chains: Sequence
    ):
        self.root_dir = DirHandle(Path(sf_dir, "pisa"), make=False)
        self.nprocs = nprocs

        self.parameters_handle = FileHandle(self.root_dir / "pisa.params")
        self.bin_path = FileHandle(self.root_dir / "pisaEnergy_linux")

        self.target_chains = tuple(target_chains)
        self.binder_chains = tuple(binder_chains)

    def __parse_output__(self, *, score_stdout=None, score_file=None) -> float:
        try:
            pisa_score = float(score_stdout)
        except ValueError as e:
            raise ValueError(f"{self} couldn't parse {score_stdout}.") from e

        return pisa_score

    def __pisa_worker__(self, frames_path: Path, i: int) -> Tuple[int, float]:
        pdb_frame = Path(frames_path, f"complex-{i}.pdb")
        # TODO: if chains are to be renamed A and B always when frames are
        # extracted (at initialize_scoring_folder()),
        # then I  should just type the 'A' and 'B' literals.
        target = "".join(self.target_chains)
        binder = "".join(self.binder_chains)

        comando_pisa = (
            f"{self.bin_path} {pdb_frame} {target} {binder} {self.parameters_handle}"
        )

        p = sp.run(comando_pisa, stdout=sp.PIPE, stderr=sp.PIPE, shell=True, text=True)
        pisa_score = self.__parse_output__(score_stdout=p.stdout)

        return i, pisa_score

    def __call__(self, *, nframes: int, frames_path: Path) -> List[float]:

        self.results_dir = DirHandle(Path(frames_path, "pisa"), make=True)
        scores: List[float] = [0] * (nframes + 1)

        with cf.ProcessPoolExecutor(max_workers=self.nprocs) as exe:
            futuros: List[cf.Future] = []
            for i in range(nframes + 1):
                futuros.append(exe.submit(self.__pisa_worker__, frames_path, i))

            timeout = self.TIMEOUT_PER_FRAME * nframes
            try:
                for futu in cf.as_completed(futuros, timeout=timeout):
                    if futu.exception():
                        logging.error(
                            f"Exception while running pisa: {futu.exception()}"
                        )
                        raise futu.exception()  # type: ignore
                    j, score = futu.result()
                    scores[j] = score
            except cf.TimeoutError as e:
                logging.error("pisa subprocess timed out.")
                raise e
        return scores
