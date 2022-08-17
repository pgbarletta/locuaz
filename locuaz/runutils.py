from attrs import define, field
from biobb_md.gromacs.grompp import grompp
from biobb_analysis.gromacs.gmx_trjconv_str import gmx_trjconv_str
from biobb_md.gromacs.mdrun import mdrun
from molecules import AbstractComplex, ZipTopology
from typing import Dict, Tuple
from pathlib import Path
from fileutils import DirHandle, FileHandle
from projectutils import WorkProject


@define
class MDrun:
    dir: DirHandle = field(converter=DirHandle)  # type: ignore
    gmx_path: str = field(converter=str, kw_only=True)
    mdp: FileHandle = field(converter=FileHandle, kw_only=True)  # type: ignore
    gpu_id: int = field(converter=int, kw_only=True, default=0)
    num_threads_omp: int = field(converter=int, kw_only=True)
    num_threads_mpi: int = field(converter=int, kw_only=True)
    dev: str = field(converter=str, kw_only=True)
    out_name: str = field(converter=str, kw_only=False)

    tpr: FileHandle = field(converter=FileHandle, init=False)  # type: ignore

    @classmethod
    def min(cls, root_dir: Path, *, work_pjct: WorkProject, out_name="min") -> "MDrun":
        obj = cls(
            root_dir,
            gmx_path=work_pjct.config["md"]["gmx_bin"],
            mdp=Path(work_pjct.mdps["min_mdp"]),
            num_threads_omp=work_pjct.config["md"]["omp_procs"],
            num_threads_mpi=work_pjct.config["md"]["mpi_procs"],
            dev="-nb gpu -pin on -pinoffset 0 -pinstride 1",
            out_name=out_name,
        )

        return obj

    @classmethod
    def nvt(
        cls, root_dir: Path, *, work_pjct: WorkProject, gpu_id: int = 0, out_name="nvt"
    ) -> "MDrun":
        obj = cls(
            root_dir,
            gmx_path=work_pjct.config["md"]["gmx_bin"],
            mdp=Path(work_pjct.mdps["nvt_mdp"]),
            gpu_id=gpu_id,
            num_threads_omp=work_pjct.config["md"]["omp_procs"],
            num_threads_mpi=work_pjct.config["md"]["mpi_procs"],
            dev="-nb gpu -pme gpu -bonded gpu -pin on -pinoffset 0 -pinstride 1 -pmefft gpu",
            out_name=out_name,
        )

        return obj

    @classmethod
    def npt(
        cls, root_dir: Path, *, work_pjct: WorkProject, gpu_id: int = 0, out_name="npt"
    ) -> "MDrun":
        obj = cls(
            root_dir,
            gmx_path=work_pjct.config["md"]["gmx_bin"],
            mdp=Path(work_pjct.mdps["npt_mdp"]),
            gpu_id=gpu_id,
            num_threads_omp=work_pjct.config["md"]["omp_procs"],
            num_threads_mpi=work_pjct.config["md"]["mpi_procs"],
            dev="-nb gpu -pme gpu -bonded gpu -pin on -pinoffset 0 -pinstride 1 -pmefft gpu",
            out_name=out_name,
        )

        return obj

    def __call__(self, complex: AbstractComplex) -> AbstractComplex:
        # Check
        self.__check_input__(complex)

        # Build .tpr file for the run
        min_tpr = Path(str(self.dir)) / (self.out_name + ".tpr")
        grompp(
            input_mdp_path=str(self.mdp),
            input_gro_path=str(complex.gro.file.path),
            input_top_zip_path=str(complex.top.file.path),
            output_tpr_path=str(min_tpr),
            properties={"gmx_path": str(self.gmx_path)},
        )

        # Run
        min_trr = Path(self.dir) / (self.out_name + ".trr")
        min_xtc = Path(self.dir) / (self.out_name + ".xtc")
        min_gro = Path(self.dir) / (self.out_name + ".gro")
        min_edr = Path(self.dir) / (self.out_name + ".edr")
        min_log = Path(self.dir) / (self.out_name + ".log")
        props = {
            "gmx_path": str(self.gmx_path),
            "num_threads_omp": self.num_threads_omp,
            "num_threads_mpi": self.num_threads_mpi,
            "gpu_id": self.gpu_id,
            "dev": self.dev,
        }

        mdrun(
            input_tpr_path=str(min_tpr),
            output_trr_path=str(min_trr),
            output_xtc_path=str(min_xtc),
            output_gro_path=str(min_gro),
            output_edr_path=str(min_edr),
            output_log_path=str(min_log),
            properties=props,
        )

        # Convert output .gro to PDB.
        min_pdb = Path(self.dir) / (self.out_name + ".pdb")
        props = {"gmx_path": str(self.gmx_path), "selection": "System"}
        gmx_trjconv_str(
            input_structure_path=str(min_gro),
            input_top_path=str(min_tpr),
            output_str_path=str(min_pdb),
            properties=props,
        )

        # Build the new complex
        new_complex = type(complex).from_complex(
            name=self.out_name,
            iter_path=self.dir,
            target_chains=complex.top.target_chains,
            binder_chains=complex.top.binder_chains,
        )

        return new_complex

    def __check_input__(self, complex: AbstractComplex) -> None:

        assert complex.dir_handle == self.dir, f"Input complex directory: "
        f"{complex.dir_handle}, does not match {type(self)}'s directory: {self.dir}"

        assert isinstance(complex.top, ZipTopology), f"Topology from  input complex "
        f"should be in zip format. Current topology: {complex.top}."
