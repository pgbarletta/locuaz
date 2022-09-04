from asyncio import as_completed
import logging
import os
import re
from pathlib import Path
from collections.abc import Sequence
import subprocess as sp
import concurrent.futures as cf
from typing import Tuple, List

from fileutils import FileHandle, DirHandle, copy_to
from abstractscoringfunction import AbstractScoringFunction


class haddock(AbstractScoringFunction):
    template_scoring_inp_handle: FileHandle
    haddock_protocols_dir: DirHandle
    haddock_toppar_dir: DirHandle
    rescoring_scripts_dir: DirHandle
    TIMEOUT_PER_FRAME: int = 30

    def __init__(
        self, sf_dir, nprocs=2, *, target_chains: Sequence, binder_chains: Sequence
    ):
        self.root_dir = DirHandle(Path(sf_dir, "haddock"), make=False)
        self.nprocs = nprocs
        
        self.bin_path = FileHandle(
            Path(self.root_dir, "cns_solve_1.3/ibm-ppc64le-linux/bin/cns")
        )
        self.template_scoring_inp_handle = Path(self.root_dir, "template_scoring.inp")
        self.haddock_protocols_dir = DirHandle(
            Path(self.root_dir, "haddock2.1", "protocols"), make=False
        )
        self.haddock_toppar_dir = DirHandle(
            Path(self.root_dir, "haddock2.1", "toppar"), make=False
        )
        self.rescoring_scripts_dir = DirHandle(
            Path(sf_dir, "rescoring-scripts"), make=False
        )
        self.target_chains = tuple(target_chains)
        self.binder_chains = tuple(binder_chains)

        # Set up environment:
        os.environ["HADDOCK"] = str(Path(self.root_dir, "haddock2.1"))
        cns_environment_script = Path(self.root_dir, "cns_solve_1.3/cns_solve_env")
        sp.run("/bin/csh " + str(cns_environment_script), shell=True)
        haddock_configure_script = Path(
            self.root_dir, "haddock2.1/haddock_configure.csh"
        )
        sp.run("/bin/csh " + str(haddock_configure_script), shell=True)

    def __pdb_chain_segid__(self, file_in: Path, file_out: Path) -> None:
        # For some reason, haddock requires the chainID to be repeated
        # onto the 73th column.
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
            #(`results_dir``), so haddock will see it.
            f.write('"' + str(mod_pdb_frame.name) + '"')

        # .inp config file for haddock.
        scorin_inp_file = Path(self.results_dir, "scoring_" + str(i) + ".inp")
        with open(self.template_scoring_inp_handle, "r") as f:
            lines = f.readlines()
        with open(scorin_inp_file, "w") as f:
            # Idem., writing only the relative path to the .pdb file.
            linea_system_name = re.sub("XYZ", str(mutant_file_list.name), lines[0])
            f.write(linea_system_name)
            
            for i in range(1, len(lines)):
                f.write(lines[i])

        return scorin_inp_file

    def __parse_output__(self, *, score_stdout=None, score_file=None) -> float:
        try:
            with open(score_file, "r") as f:
                lineas = f.readlines()
                score_haddock = float(lineas[8].split()[1][0:-1])
        except ValueError as e:
            raise ValueError(f"{self} couldn't parse {score_file}.") from e

        return score_haddock

    def __haddock_worker__(self, frames_path: Path, i: int) -> Tuple[int, float]:
        
        pdb_frame = Path(frames_path, f"complex-{i}.pdb")
        mod_pdb_frame = Path(self.results_dir, f"mod_complex-{i}.pdb")
        self.__pdb_chain_segid__(pdb_frame, mod_pdb_frame)
        scorin_inp_file = self.__init_frame_scoring__(i, mod_pdb_frame)
        output = Path(self.results_dir, f"log_{i}.out")
        
        comando_haddock = f"{self.bin_path} <  {scorin_inp_file} > {output}"
        sp.run(comando_haddock,stdout=sp.PIPE,stderr=sp.PIPE,cwd=self.results_dir,shell=True,text=True)

        output_haddock_file = Path(self.results_dir, f"mod_complex-{i}_conv.psf")
        score_haddock = self.__parse_output__(score_file=output_haddock_file)        

        return i, score_haddock

    def __call__(self, *, nframes: int, frames_path: Path) -> List[float]:

        self.results_dir = DirHandle(Path(frames_path, "haddock"), make=True)
        self.__initialize_scoring_dir__()
        scores: List[float] = [0] * (nframes + 1)
        with cf.ProcessPoolExecutor(max_workers=self.nprocs) as exe:
            futuros: List[cf.Future] = []
            for i in range(nframes+1):
                futuros.append(exe.submit(self.__haddock_worker__, frames_path, i))

            timeout = self.TIMEOUT_PER_FRAME * nframes
            try:
                for futu in cf.as_completed(futuros, timeout=timeout):
                    if futu.exception():
                        logging.error(
                            f"Exception while running hadock: {futu.exception()}"
                        )
                        raise futu.exception()  # type: ignore
                    j, score = futu.result()
                    scores[j] = score
            except cf.TimeoutError as e:
                logging.error("haddock subprocess timed out.")
                raise e

        return scores


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

    