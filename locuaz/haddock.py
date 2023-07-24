import concurrent.futures as cf
import os
import re
import subprocess as sp
from pathlib import Path
from typing import Tuple, List, Union

from locuaz.abstractscorer import AbstractScorer
from locuaz.complex import GROComplex
from locuaz.fileutils import FileHandle, DirHandle, copy_to


class Haddock(AbstractScorer):
    template_scoring_inp_handle: FileHandle
    haddock_protocols_dir: DirHandle
    haddock_toppar_dir: DirHandle
    rescoring_scripts_dir: DirHandle
    TIMEOUT_PER_FRAME: int = 30

    def __init__(self, sf_dir, *, nthreads=2, mpi_procs=2) -> None:
        super().__init__(sf_dir, nthreads=nthreads, mpi_procs=mpi_procs)
        self.template_scoring_inp_handle = FileHandle(
            Path(self.root_dir, "template_scoring.inp")
        )
        self.haddock_protocols_dir = DirHandle(
            Path(self.root_dir, "protocols"), make=False
        )
        self.haddock_toppar_dir = DirHandle(Path(self.root_dir, "toppar"), make=False)
        self.rescoring_scripts_dir = DirHandle(
            Path(self.root_dir, "rescoring-scripts"), make=False
        )

        # Set up environment:
        os.environ["HADDOCK"] = str(Path(self.root_dir))

        cns_environment_script = Path(self.root_dir, "cns_solve_env")
        sp.run("/bin/csh " + str(cns_environment_script), shell=True)

        haddock_configure_script = Path(self.root_dir, "haddock_configure.csh")
        sp.run("/bin/csh " + str(haddock_configure_script), shell=True)

    @staticmethod
    def __pdb_chain_segid__(file_in: Path, file_out: Path) -> None:
        # For some reason, haddock requires the chainID to be repeated
        # onto the 73rd column as a segid.
        with open(file_in, "r") as sources:
            lines = sources.readlines()
        with open(file_out, "w") as sources:
            for linea in lines:
                if linea[0:4] == "ATOM":
                    line = linea[0:72] + linea[21] + linea[73:]
                else:
                    line = linea
                sources.write(line)

    def __init_frame_scoring__(self, i: int, mod_pdb_frame: Path) -> Path:
        # .list file with the path to the PDB
        mutant_file_list = Path(self.results_dir, "pdb_to_haddock_" + str(i) + ".list")
        with open(mutant_file_list, "w") as f:
            # Writing only the name of `mod_pdb_frame` because haddock can't deal
            # with long filenames. The subprocess will be run from the same cwd dir
            # (`results_dir``), so haddock will see it.
            f.write('"' + str(mod_pdb_frame.name) + '"')

        # .inp config file for haddock.
        scoring_inp_file = Path(self.results_dir, "scoring_" + str(i) + ".inp")
        with open(self.template_scoring_inp_handle, "r") as f:
            lines = f.readlines()
        with open(scoring_inp_file, "w") as f:
            # Idem., writing only the relative path to the .pdb file.
            linea_system_name = re.sub("XYZ", str(mutant_file_list.name), lines[0])
            f.write(linea_system_name)

            for i in range(1, len(lines)):
                f.write(lines[i])

        return scoring_inp_file

    def __parse_outfile_(self, score_file: Union[Path, FileHandle], original_command: str) -> float:
        assert (
            score_file is not None
        ), f"This shouldn't happen. {self} couldn't parse {score_file}\nfrom: \n{original_command}"
        try:
            with open(score_file, "r") as f:
                lineas = f.readlines()
                score_haddock = float(lineas[8].split()[1][0:-1])
        except (ValueError, IndexError) as e:
            raise ValueError(
                f"{self} couldn't parse {score_file}\nfrom: \n{original_command}"
            ) from e

        return score_haddock

    def __haddock_worker__(self, frames_path: Path, i: int) -> Tuple[int, float]:

        pdb_frame = Path(frames_path, f"complex-{i}.pdb")
        mod_pdb_frame = Path(self.results_dir, f"mod_complex-{i}.pdb")
        self.__pdb_chain_segid__(pdb_frame, mod_pdb_frame)
        scorin_inp_file = self.__init_frame_scoring__(i, mod_pdb_frame)
        output = Path(self.results_dir, f"log_{i}.out")

        comando_haddock = f"{self.bin_path} <  {scorin_inp_file} > {output}"
        p = sp.run(
            comando_haddock,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            cwd=self.results_dir,
            shell=True,
            text=True,
        )

        output_haddock_file = Path(self.results_dir, f"mod_complex-{i}_conv.psf")
        self.__assert_scorer_outfile__(output_haddock_file, stdout=p.stdout, stderr=p.stderr,
                                                 command=comando_haddock)
        score_haddock = self.__parse_outfile_(output_haddock_file, comando_haddock)

        return i, score_haddock

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
        self.__initialize_scoring_dir__()
        # The first unused frames will be discarded later
        scores: List[float] = [0] * end
        with cf.ProcessPoolExecutor(max_workers=self.nthreads) as exe:
            futuros: List[cf.Future] = []
            for i in range(start, end):
                futuros.append(exe.submit(self.__haddock_worker__, frames_path, i))

            timeout = self.TIMEOUT_PER_FRAME * nframes
            try:
                for futu in cf.as_completed(futuros, timeout=timeout):
                    if futu.exception():
                        print(
                            f"Exception while running hadock: {futu.exception()}",
                            flush=True,
                        )
                        raise futu.exception()  # type: ignore
                    j, score = futu.result()
                    scores[j] = score
            except cf.TimeoutError as e:
                print(
                    f"{self.name} subprocess timed out after {timeout} seconds.",
                    flush=True,
                )
                raise e

        # Discard the first 0 frames
        return scores[start:]

    def __initialize_scoring_dir__(self) -> None:

        copy_to(
            FileHandle(Path(self.rescoring_scripts_dir, "ligand.param")),
            self.results_dir,
        )
        copy_to(
            FileHandle(Path(self.rescoring_scripts_dir, "ligand.top")), self.results_dir
        )
        copy_to(
            FileHandle(Path(self.rescoring_scripts_dir, "calc_free-ene.cns")),
            self.results_dir,
        )
        copy_to(
            FileHandle(Path(self.rescoring_scripts_dir, "def_solv_param.cns")),
            self.results_dir,
        )
        copy_to(
            FileHandle(Path(self.rescoring_scripts_dir, "flex_segment_back.cns")),
            self.results_dir,
        )
        copy_to(
            FileHandle(Path(self.rescoring_scripts_dir, "print_coorheader.cns")),
            self.results_dir,
        )
        copy_to(
            FileHandle(Path(self.rescoring_scripts_dir, "scale_inter_only.cns")),
            self.results_dir,
        )
        copy_to(
            FileHandle(Path(self.rescoring_scripts_dir, "scale_intra_only.cns")),
            self.results_dir,
        )
        copy_to(
            FileHandle(Path(self.rescoring_scripts_dir, "separate.cns")),
            self.results_dir,
        )
        copy_to(
            FileHandle(Path(self.rescoring_scripts_dir, "symmultimer.cns")),
            self.results_dir,
        )
