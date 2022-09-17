import logging
from textwrap import wrap
from typing import Optional
from attrs import define, field
from pathlib import Path
from shutil import SameFileError

from biobb_md.gromacs.grompp import Grompp
from biobb_analysis.gromacs.gmx_trjconv_str import GMXTrjConvStr
from biobb_md.gromacs.mdrun import Mdrun
from biobb_analysis.gromacs.gmx_image import GMXImage
from molecules import AbstractComplex, ZipTopology, copy_mol_to

from fileutils import DirHandle, FileHandle, copy_to
from projectutils import WorkProject
from primitives import launch_biobb


@define(frozen=True)
class MDrun:
    dir: DirHandle = field(converter=DirHandle)  # type: ignore
    gmx_path: str = field(converter=str, kw_only=True)
    mdp: FileHandle = field(converter=FileHandle, kw_only=True)  # type: ignore
    gpu_id: int = field(converter=int, kw_only=True, default=0)
    pinoffset: int = field(converter=int, kw_only=True, default=0)
    num_threads_omp: int = field(converter=int, kw_only=True)
    num_threads_mpi: int = field(converter=int, kw_only=True)
    dev: str = field(converter=str, kw_only=True)
    out_name: str = field(converter=str, kw_only=True)
    nojump: bool = field(converter=bool, kw_only=True, default=False)

    @classmethod
    def min(
        cls,
        root_dir: Path,
        *,
        work_pjct: WorkProject,
        gpu_id: int = 0,
        pinoffset: int = 0,
        out_name="min",
    ) -> "MDrun":

        obj = cls(
            root_dir,
            gmx_path=work_pjct.config["md"]["gmx_bin"],
            mdp=Path(work_pjct.mdps["min_mdp"]),
            gpu_id=gpu_id,
            pinoffset=pinoffset,
            num_threads_omp=work_pjct.config["md"]["omp_procs"],
            num_threads_mpi=work_pjct.config["md"]["mpi_procs"],
            dev=f"-nb gpu -pin on -pinoffset {pinoffset}",
            out_name=out_name,
        )

        return obj

    @classmethod
    def nvt(
        cls,
        root_dir: Path,
        *,
        work_pjct: WorkProject,
        gpu_id: int = 0,
        pinoffset: int = 0,
        out_name="nvt",
    ) -> "MDrun":

        obj = cls(
            root_dir,
            gmx_path=work_pjct.config["md"]["gmx_bin"],
            mdp=Path(work_pjct.mdps["nvt_mdp"]),
            gpu_id=gpu_id,
            pinoffset=pinoffset,
            num_threads_omp=work_pjct.config["md"]["omp_procs"],
            num_threads_mpi=work_pjct.config["md"]["mpi_procs"],
            dev=f"-nb gpu -pme gpu -bonded gpu -pmefft gpu -pin on -pinoffset {pinoffset}",
            out_name=out_name,
        )

        return obj

    @classmethod
    def npt(
        cls,
        root_dir: Path,
        *,
        work_pjct: WorkProject,
        gpu_id: int = 0,
        pinoffset: int = 0,
        out_name="npt",
        nojump=True,
    ) -> "MDrun":

        obj = cls(
            root_dir,
            gmx_path=work_pjct.config["md"]["gmx_bin"],
            mdp=Path(work_pjct.mdps["npt_mdp"]),
            gpu_id=gpu_id,
            pinoffset=pinoffset,
            num_threads_omp=work_pjct.config["md"]["omp_procs"],
            num_threads_mpi=work_pjct.config["md"]["mpi_procs"],
            dev=f"-nb gpu -pme gpu -bonded gpu -pmefft gpu -pin on -pinoffset {pinoffset}",
            out_name=out_name,
            nojump=nojump,
        )

        return obj

    def __call__(self, complex: AbstractComplex) -> AbstractComplex:
        # Check
        self.__check_input__(complex)

        # Build .tpr file for the run
        run_tpr = Path(str(self.dir)) / (self.out_name + ".tpr")
        grompepe = Grompp(
            input_mdp_path=str(self.mdp),
            input_gro_path=str(complex.gro.file.path),
            input_top_zip_path=str(complex.top.file.path),
            input_cpt_path=str(complex.cpt),
            output_tpr_path=str(run_tpr),
            properties={"gmx_path": str(self.gmx_path)},
        )
        launch_biobb(grompepe)

        # Run
        run_trr = Path(self.dir) / (self.out_name + ".trr")
        run_xtc = Path(self.dir) / (self.out_name + ".xtc")
        run_gro = Path(self.dir) / (self.out_name + ".gro")
        run_edr = Path(self.dir) / (self.out_name + ".edr")
        run_log = Path(self.dir) / (self.out_name + ".log")
        run_cpt = Path(self.dir) / (self.out_name + ".cpt")
        props = {
            "gmx_path": str(self.gmx_path),
            "num_threads_omp": self.num_threads_omp,
            "num_threads_mpi": self.num_threads_mpi,
            "gpu_id": self.gpu_id,
            "dev": self.dev,
        }

        runner = Mdrun(
            input_tpr_path=str(run_tpr),
            input_cpt_path=str(complex.cpt),
            output_trr_path=str(run_trr),
            output_xtc_path=str(run_xtc),
            output_gro_path=str(run_gro),
            output_edr_path=str(run_edr),
            output_log_path=str(run_log),
            output_cpt_path=str(run_cpt),
            properties=props,
        )
        launch_biobb(runner)

        if self.nojump:
            tmp_xtc = Path(self.dir) / "tmp.xtc"

            wrap_mol = GMXImage(
                input_traj_path=str(run_xtc),
                input_top_path=str(complex.gro),
                output_traj_path=str(tmp_xtc),
                properties={
                    "gmx_path": str(self.gmx_path),
                    "center_selection": "Protein",
                    "output_selection": "System",
                    "pbc": "nojump",
                },
            )
            launch_biobb(wrap_mol)

            copy_to(FileHandle(tmp_xtc), Path(self.dir), name=run_xtc.name)
            tmp_xtc.unlink()

        # Finally, build the Complex.
        try:
            copy_mol_to(complex.top, self.dir, self.out_name + ".zip")
        except SameFileError:
            logging.warning(
                f"Attempted to run MD on a finished run, starting from: {complex.cpt} "
            )

        new_complex = type(complex).from_gro_zip(
            name=self.out_name,
            input_dir=self.dir,
            target_chains=complex.top.target_chains,
            binder_chains=complex.top.binder_chains,
            gmx_bin=self.gmx_path,
        )

        # Convert output .gro to PDB.
        run_pdb = Path(self.dir) / (self.out_name + ".pdb")
        props = {"gmx_path": str(self.gmx_path), "selection": "System"}
        gro_to_pdb = GMXTrjConvStr(
            input_structure_path=str(run_gro),
            input_top_path=str(run_tpr),
            output_str_path=str(run_pdb),
            properties=props,
        )
        launch_biobb(gro_to_pdb)

        # Build the new complex
        new_complex = type(complex).from_complex(
            name=self.out_name,
            iter_path=self.dir,
            target_chains=complex.top.target_chains,
            binder_chains=complex.top.binder_chains,
            gmx_bin=self.gmx_path,
        )

        return new_complex

    def __check_input__(self, complex: AbstractComplex) -> None:

        assert complex.dir == self.dir, f"Input complex directory: "
        f"{complex.dir}, does not match {type(self)}'s directory: {self.dir}"

        assert isinstance(complex.top, ZipTopology), f"Topology from  input complex "
        f"should be in zip format. Current topology: {complex.top}."
