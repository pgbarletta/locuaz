from pathlib import Path
from shutil import SameFileError
from typing import Tuple
import subprocess as sp

from attrs import define, field
from biobb_analysis.gromacs.gmx_trjconv_str import GMXTrjConvStr
from biobb_gromacs.gromacs.grompp import Grompp
from biobb_gromacs.gromacs.mdrun import Mdrun

from complex import AbstractComplex, GROComplex
from fileutils import DirHandle, FileHandle
from molecules import ZipTopology, copy_mol_to
from fixbox import fix_box_cpx
from primitives import launch_biobb, GromacsError
from projectutils import WorkProject


@define(frozen=True)
class MDrun:
    dir: DirHandle = field(converter=DirHandle)  # type: ignore
    binary_path: str = field(converter=str, kw_only=True)
    mdp: FileHandle = field(converter=FileHandle, kw_only=True)  # type: ignore
    gpu_id: int = field(converter=int, kw_only=True, default=0)
    pinoffset: int = field(converter=int, kw_only=True, default=0)
    num_threads_omp: int = field(converter=int, kw_only=True)
    num_threads_mpi: int = field(converter=int, kw_only=True)
    dev: str = field(converter=str, kw_only=True)
    out_name: str = field(converter=str, kw_only=True)
    image_after: bool = field(converter=bool, kw_only=True, default=False)

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
            DirHandle(root_dir),
            binary_path=work_pjct.config["md"]["gmx_bin"],
            mdp=FileHandle(Path(work_pjct.mdps["min_mdp"])),
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
            DirHandle(root_dir),
            binary_path=work_pjct.config["md"]["gmx_bin"],
            mdp=FileHandle(Path(work_pjct.mdps["nvt_mdp"])),
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
        image_after=True,
    ) -> "MDrun":

        obj = cls(
            DirHandle(root_dir),
            binary_path=work_pjct.config["md"]["gmx_bin"],
            mdp=FileHandle(Path(work_pjct.mdps["npt_mdp"])),
            gpu_id=gpu_id,
            pinoffset=pinoffset,
            num_threads_omp=work_pjct.config["md"]["omp_procs"],
            num_threads_mpi=work_pjct.config["md"]["mpi_procs"],
            dev=f"-nb gpu -pme gpu -bonded gpu -pmefft gpu -pin on -pinoffset {pinoffset}",
            out_name=out_name,
            image_after=image_after,
        )

        return obj

    def __call__(self, complex: GROComplex) -> Tuple[bool, GROComplex]:
        # Check
        self.__check_input__(complex)

        # Build .tpr file for the run
        run_tpr = Path(str(self.dir)) / (self.out_name + ".tpr")
        if complex.cpt:
            grompepe = Grompp(
                input_mdp_path=str(self.mdp),
                input_gro_path=str(complex.gro.file.path),
                input_top_zip_path=str(complex.top.file.path),
                input_cpt_path=str(complex.cpt),
                input_ndx_path=str(complex.ndx),
                output_tpr_path=str(run_tpr),
                properties={"binary_path": str(self.binary_path)},
            )
        else:
            grompepe = Grompp(
                input_mdp_path=str(self.mdp),
                input_gro_path=str(complex.gro.file.path),
                input_top_zip_path=str(complex.top.file.path),
                input_ndx_path=str(complex.ndx),
                output_tpr_path=str(run_tpr),
                properties={"binary_path": str(self.binary_path)},
            )
        launch_biobb(grompepe, backup_dict=Path(self.dir))

        # Run
        run_trr = Path(self.dir) / (self.out_name + ".trr")
        run_xtc = Path(self.dir) / (self.out_name + ".xtc")
        run_gro = Path(self.dir) / (self.out_name + ".gro")
        run_edr = Path(self.dir) / (self.out_name + ".edr")
        run_log = Path(self.dir) / (self.out_name + ".log")
        run_cpt = Path(self.dir) / (self.out_name + ".cpt")
        run_pux = Path(self.dir) / (f"pullx_{self.out_name}.xvg")
        run_puf = Path(self.dir) / (f"pullf_{self.out_name}.xvg")
        props = {
            "binary_path": str(self.binary_path),
            "num_threads_omp": self.num_threads_omp,
            "num_threads_mpi": self.num_threads_mpi,
            "gpu_id": self.gpu_id,
            "dev": f"{self.dev} -px {run_pux} -pf {run_puf}",
        }

        comando_md = f'{self.binary_path} -nobackup -nocopyright mdrun'
        comando_md += f' -s {run_tpr} -c {run_gro}'
        comando_md += f' -cpi {run_cpt}' if run_cpt.is_file() else ""
        comando_md += f' -px {run_pux} -pf {run_puf}'
        comando_md += f' -o {run_trr} -x {run_xtc} -g {run_log} -cpo {run_cpt}'
        comando_md += f' {self.dev} -gpu_id {self.gpu_id} -ntmpi {self.num_threads_mpi} -ntomp {self.num_threads_omp}'

        try:
            p = sp.run(
                comando_md,
                stdout=sp.PIPE,
                stderr=sp.PIPE,
                cwd=Path(self.dir),
                shell=True,
                text=True,
            )
        except (RuntimeError, Exception):
            raise GromacsError(f"MDRun error. stdout:\n{p.stdout}\n"f"stderr:\n{p.stderr}")

        # Finally, build the Complex.
        try:
            copy_mol_to(complex.top, self.dir, f"{self.out_name}.zip")
        except SameFileError:
            # This happens when there was an attempt to run MD on a finished run
            # and gromacs started from `name`.cpt. This shouldn't happen anymore.
            print(f"{complex.top} same as {self.out_name}.zip", flush=True)
            pass

        new_complex = type(complex).from_gro_zip(
            name=self.out_name,
            input_dir=Path(self.dir),
            target_chains=complex.top.target_chains,
            binder_chains=complex.top.binder_chains,
            gmx_bin=self.binary_path,
        )

        all_atoms_in_box = True
        if self.image_after:
            run_pdb = Path(self.dir, f"{self.out_name}.pdb")
            all_atoms_in_box, _ = fix_box_cpx(new_complex, run_pdb, str(self.binary_path))
        else:
            # Convert output .gro to PDB.
            run_pdb = Path(self.dir) / (self.out_name + ".pdb")
            props = {"binary_path": str(self.binary_path), "selection": "System"}
            gro_to_pdb = GMXTrjConvStr(
                input_structure_path=str(run_gro),
                input_top_path=str(run_tpr),
                output_str_path=str(run_pdb),
                properties=props,
            )
            launch_biobb(gro_to_pdb, backup_dict=Path(self.dir))

        # Build the new complex
        new_complex = type(complex).from_complex(
            name=self.out_name,
            iter_path=Path(self.dir),
            target_chains=complex.top.target_chains,
            binder_chains=complex.top.binder_chains,
            gmx_bin=self.binary_path,
        )

        return all_atoms_in_box, new_complex

    def __check_input__(self, complex: AbstractComplex) -> None:

        assert complex.dir == self.dir, f"Input complex directory: "
        f"{complex.dir}, does not match {type(self)}'s directory: {self.dir}"

        assert isinstance(complex.top, ZipTopology), f"Topology from input complex "
        f"should be in zip format. Current topology: {complex.top}."
