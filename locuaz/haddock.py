from pathlib import Path
from collections import Sequence
import subprocess as sp
from fileutils import FileHandle, DirHandle
from abstractscoringfunction import AbstractScoringFunction
from collections import Sequence
from shutil import copyfile
import os
import re


class haddock(AbstractScoringFunction):
    CPU_TO_MEM_RATIO: int = 8
    template_scoring_inp_handle: FileHandle
    haddock_protocols_dir: DirHandle
    haddock_toppar_dir: DirHandle
    rescoring_scripts_dir: DirHandle

    def __init__(
        self, sf_dir, nprocs=2, *, target_chains: Sequence, binder_chains: Sequence
    ):
        self.root_dir = DirHandle(Path(sf_dir, "haddock"), make=False)
        self.nprocs = nprocs
        self.max_concurrent_jobs = self.CPU_TO_MEM_RATIO * self.nprocs

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
            f.write('"' + str(mod_pdb_frame) + '"')

        # .inp config file for haddock.
        scorin_inp_file = Path(self.results_dir, "scoring_" + str(i) + ".inp")
        with open(self.template_scoring_inp_handle, "r") as f:
            lines = f.readlines()
        with open(scorin_inp_file, "w") as f:
            linea_system_name = re.sub("XYZ", str(mutant_file_list), lines[0])

            ncomponents = len(self.target_chains) + len(self.binder_chains)
            linea_ncomponents = re.sub(
                "ncomponents=2", f"ncomponents={ncomponents}", lines[1]
            )
            # TODO: check con Miguel
            segid_1 = "".join(self.target_chains)
            linea_prot_segid_1 = re.sub(
                'prot_segid_1="A"', f'prot_segid_1="{segid_1}"', lines[2]
            )
            segid_2 = "".join(self.binder_chains)
            linea_prot_segid_2 = re.sub(
                'prot_segid_2="B"', f'prot_segid_2="{segid_2}"', lines[3]
            )
            f.write(linea_system_name)
            f.write(linea_ncomponents)
            f.write(linea_prot_segid_1)
            f.write(linea_prot_segid_2)
            for i in range(4, len(lines)):
                f.write(lines[i])

        return scorin_inp_file

    def __submit_batch__(
        self,
        start: int,
        stop: int,
        frames_path: Path,
    ):
        processos = []
        for i in range(start, stop):
            pdb_frame = Path(frames_path, "complex-" + str(i) + ".pdb")
            mod_pdb_frame = Path(self.results_dir, "mod_complex-" + str(i) + ".pdb")
            self.__pdb_chain_segid__(pdb_frame, mod_pdb_frame)
            scorin_inp_file = self.__init_frame_scoring__(i, mod_pdb_frame)

            comando_haddock = (
                str(self.bin_path)
                + " < "
                + str(scorin_inp_file)
                + " > "
                + str(Path(self.results_dir, "log_" + str(i) + ".out"))
            )

            processos.append(
                sp.Popen(
                    comando_haddock,
                    stdout=sp.PIPE,
                    stderr=sp.PIPE,
                    cwd=self.results_dir,
                    shell=True,
                    text=True,
                )
            )

        scores_haddock = []
        for i, proc in enumerate(processos):
            _, _ = proc.communicate()
            output_haddock_file = Path(
                self.results_dir, "mod_complex-" + str(i) + "_conv.psf"
            )
            with open(output_haddock_file, "r") as f:
                lineas = f.readlines()
                scores_haddock.append(float(lineas[8].split()[1][0:-1]))

        return scores_haddock

    def __initialize_scoring_dir__(self) -> None:
        copyfile(
            Path(self.rescoring_scripts_dir, "ligand.param"),
            Path(self.results_dir, "ligand.param"),
        )
        copyfile(
            Path(self.rescoring_scripts_dir, "ligand.top"),
            Path(self.results_dir, "ligand.top"),
        )

        copyfile(
            Path(self.rescoring_scripts_dir, "calc_free-ene.cns"),
            Path(self.results_dir, "calc_free-ene.cns"),
        )
        copyfile(
            Path(self.rescoring_scripts_dir, "def_solv_param.cns"),
            Path(self.results_dir, "def_solv_param.cns"),
        )
        copyfile(
            Path(self.rescoring_scripts_dir, "flex_segment_back.cns"),
            Path(self.results_dir, "flex_segment_back.cns"),
        )
        copyfile(
            Path(self.rescoring_scripts_dir, "print_coorheader.cns"),
            Path(self.results_dir, "print_coorheader.cns"),
        )
        copyfile(
            Path(self.rescoring_scripts_dir, "scale_inter_only.cns"),
            Path(self.results_dir, "scale_inter_only.cns"),
        )
        copyfile(
            Path(self.rescoring_scripts_dir, "scale_intra_only.cns"),
            Path(self.results_dir, "scale_intra_only.cns"),
        )
        copyfile(
            Path(self.rescoring_scripts_dir, "separate.cns"),
            Path(self.results_dir, "separate.cns"),
        )
        copyfile(
            Path(self.rescoring_scripts_dir, "symmultimer.cns"),
            Path(self.results_dir, "symmultimer.cns"),
        )

    def __call__(self, *, nframes: int, frames_path: Path):
        print(" -- HADDOCK scoring -- ")

        self.results_dir = DirHandle(Path(frames_path, "haddock"), make=True)
        self.__initialize_scoring_dir__()
        steps = list(range(0, nframes + 1, self.max_concurrent_jobs))
        scores = []
        for start, stop in zip(steps[0::1], steps[1::1]):
            step_scores = self.__submit_batch__(start, stop, frames_path)
            scores += step_scores

        # Remaining frames:
        step_scores = self.__submit_batch__(steps[-1], nframes + 1, frames_path)
        scores += step_scores

        return scores
