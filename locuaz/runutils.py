import os
from pathlib import Path
import shutil
from typing import Tuple, Union, Dict, List, Final, Optional
import subprocess as sp
from math import ceil
import numpy as np
from itertools import chain
from zipfile import ZipFile
from warnings import warn

from attrs import define, field
from biobb_analysis.gromacs.gmx_trjconv_str import GMXTrjConvStr
from biobb_gromacs.gromacs.grompp import Grompp

from locuaz.complex import AbstractComplex, GROComplex
from locuaz.fileutils import DirHandle, FileHandle
from locuaz.molecules import ZipTopology, copy_mol_to
from locuaz.moleculesutils import set_posres
from locuaz.fixbox import fix_box_cpx
from locuaz.primitives import launch_biobb, GromacsError
from locuaz.projectutils import Epoch

MAX_OMP_THREADS: Final[int] = 16
OMP_THREADS: Final[List[int]] = list(range(2, MAX_OMP_THREADS + 1))
OMP_MPI_THREADS: Final[List[int]] = list(range(4, MAX_OMP_THREADS + 1 + 2))


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
    maxwarn: int = field(converter=int, kw_only=True, default=0)
    restraints: Optional[Dict[str, float]] = field(kw_only=True, default=None)

    @classmethod
    def min(
        cls,
        root_dir: Union[Path, DirHandle],
        *,
        gmx_mdrun: str,
        min_mdp: Path,
        gpu_id: int = 0,
        omp_threads: int = 2,
        mpi_threads: int = 1,
        pinoffset: int = 0,
        out_name="min",
        maxwarn: int = 0,
    ) -> "MDrun":
        obj = cls(
            DirHandle(root_dir),
            binary_path=gmx_mdrun,
            mdp=FileHandle(min_mdp),
            gpu_id=gpu_id,
            num_threads_omp=omp_threads,
            num_threads_mpi=mpi_threads,
            pinoffset=pinoffset,
            dev=f"-nb gpu -pin on -pinoffset {pinoffset}",
            out_name=out_name,
            maxwarn=maxwarn,
        )
        return obj

    @classmethod
    def nvt(
        cls,
        root_dir: Union[Path, DirHandle],
        *,
        gmx_mdrun: str,
        nvt_mdp: Path,
        gpu_id: int = 0,
        omp_threads: int = 2,
        mpi_threads: int = 1,
        pinoffset: int = 0,
        out_name="nvt",
        maxwarn: int = 0,
    ) -> "MDrun":
        obj = cls(
            DirHandle(root_dir),
            binary_path=gmx_mdrun,
            mdp=FileHandle(nvt_mdp),
            gpu_id=gpu_id,
            num_threads_omp=omp_threads,
            num_threads_mpi=mpi_threads,
            pinoffset=pinoffset,
            dev=f"-nb gpu -pme gpu -pin on -pinoffset {pinoffset}",
            out_name=out_name,
            maxwarn=maxwarn,
        )
        return obj

    @classmethod
    def npt(
        cls,
        root_dir: Union[Path, DirHandle],
        *,
        gmx_mdrun: str,
        npt_mdp: Path,
        gpu_id: int = 0,
        omp_threads: int = 2,
        mpi_threads: int = 1,
        pinoffset: int = 0,
        out_name="npt",
        image_after=True,
        maxwarn: int = 0,
        restraints: Optional[Dict[str, float]] = None,
    ) -> "MDrun":
        obj = cls(
            DirHandle(root_dir),
            binary_path=gmx_mdrun,
            mdp=FileHandle(npt_mdp),
            gpu_id=gpu_id,
            num_threads_omp=omp_threads,
            num_threads_mpi=mpi_threads,
            pinoffset=pinoffset,
            dev=f"-nb gpu -pme gpu -pin on -pinoffset {pinoffset}",
            out_name=out_name,
            maxwarn=maxwarn,
            image_after=image_after,
            restraints=restraints,
        )
        return obj

    def __call__(self, complex: GROComplex) -> Tuple[bool, GROComplex]:
        self.__check_input__(complex)

        restrained_top = self.__set_restraints__(complex.top)

        # Build .tpr file for the run
        run_tpr = Path(str(self.dir)) / (self.out_name + ".tpr")
        if complex.cpt:
            grompepe = Grompp(
                input_mdp_path=str(self.mdp),
                input_gro_path=str(complex.gro.file.path),
                input_top_zip_path=str(restrained_top.file.path),
                input_cpt_path=str(complex.cpt),
                input_ndx_path=str(complex.ndx),
                output_tpr_path=str(run_tpr),
                properties={"binary_path": "gmx", "maxwarn": self.maxwarn},
            )
        else:
            grompepe = Grompp(
                input_mdp_path=str(self.mdp),
                input_gro_path=str(complex.gro.file.path),
                input_top_zip_path=str(restrained_top.file.path),
                input_ndx_path=str(complex.ndx),
                output_tpr_path=str(run_tpr),
                properties={"binary_path": "gmx", "maxwarn": self.maxwarn},
            )
        launch_biobb(grompepe, backup_dict=Path(self.dir))

        # Run
        run_trr = Path(self.dir) / (self.out_name + ".trr")
        run_xtc = Path(self.dir) / (self.out_name + ".xtc")
        run_gro = Path(self.dir) / (self.out_name + ".gro")
        run_edr = Path(self.dir) / (self.out_name + ".edr")
        run_log = Path(self.dir) / (self.out_name + ".log")
        run_cpt = Path(self.dir) / (self.out_name + ".cpt")
        run_pux = Path(self.dir) / f"pullx_{self.out_name}.xvg"
        run_puf = Path(self.dir) / f"pullf_{self.out_name}.xvg"

        comando_md = f"{self.binary_path} -nobackup -nocopyright"
        comando_md += f" -s {run_tpr} -c {run_gro}"
        comando_md += f" -cpi {run_cpt}" if run_cpt.is_file() else ""
        comando_md += f" -px {run_pux} -pf {run_puf} -e {run_edr}"
        comando_md += f" -o {run_trr} -x {run_xtc} -g {run_log} -cpo {run_cpt}"
        comando_md += f" {self.dev} -gpu_id {self.gpu_id}"
        comando_md += f" -ntmpi {self.num_threads_mpi} -ntomp {self.num_threads_omp}"

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
            # noinspection PyUnboundLocalVariable
            raise GromacsError(
                f"MDRun error. stdout:\n{p.stdout}\n" f"stderr:\n{p.stderr}"
            )

        self.assert_outfile(
            run_gro, stdout=p.stdout, stderr=p.stderr, command=comando_md
        )

        # Finally, build the Complex.
        try:
            copy_mol_to(complex.top, self.dir, f"{self.out_name}.zip")
        except shutil.SameFileError:
            # This happens when there was an attempt to run MD on a finished run and gromacs started from `name`.cpt.
            # This only happens when a run was interrupted after some branches were finished and others weren't.
            # Otherwise, the `npt_done` flag of the epoch would be set and this step would be skipped.
            run_dir = Path(self.dir).name
            print(
                f"Run from {run_dir} had already finished. {complex.top} same as {self.out_name}.zip",
                flush=True,
            )
            pass

        new_complex = type(complex).from_gro_zip(
            name=self.out_name,
            input_dir=Path(self.dir),
            target_chains=complex.top.target_chains,
            binder_chains=complex.top.binder_chains,
            gmx_bin="gmx",
        )
        all_atoms_in_box = True
        if self.image_after:
            run_pdb = Path(self.dir, f"{self.out_name}.pdb")
            all_atoms_in_box, _ = fix_box_cpx(new_complex, run_pdb)
        else:
            # Convert output .gro to PDB.
            run_pdb = Path(self.dir) / (self.out_name + ".pdb")
            props = {"binary_path": "gmx", "selection": "System"}
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
            branch_path=Path(self.dir),
            target_chains=complex.top.target_chains,
            binder_chains=complex.top.binder_chains,
            gmx_bin="gmx",
        )

        return all_atoms_in_box, new_complex

    def __check_input__(self, complex: AbstractComplex) -> None:
        assert complex.dir == self.dir, f"Input complex directory: "
        f"{complex.dir}, does not match {type(self)}'s directory: {self.dir}"

        assert isinstance(complex.top, ZipTopology), f"Topology from input complex "
        f"should be in zip format. Current topology: {complex.top}."

    def assert_outfile(
        self,
        out_file: Union[str, Path, FileHandle],
        *,
        stdout: str,
        stderr: str,
        command: str,
    ) -> Path:
        out_file_path = Path(out_file)
        assert out_file_path.is_file(), f"""{self} error. Can't parse: {out_file_path}
from:
{command}
with stdout:
{stdout}
and stderr:
{stderr}
"""
        return out_file_path

    def __set_restraints__(self, topology: ZipTopology) -> ZipTopology:
        if self.restraints is None:
            return topology
        zipped_top = ZipFile(topology.file.path)
        filelist = [file.filename for file in zipped_top.filelist]
        posres_filelist = {filename for filename in filelist if "posre" in filename}
        if len(posres_filelist) == 0:
            warn(
                "ZipTopology does not contain any file beginning with 'posres_...', MDrun can't set restraints."
            )
            return topology

        with zipped_top as sipesipe:
            sipesipe.extractall(Path(self.dir))

        for posre_file in posres_filelist:
            posre_strings = set(posre_file.split(".")[0].split("_"))

            if len(posre_strings.intersection(set(topology.target_chains))) != 0:
                self.__overwrite_posres_file__(
                    Path(self.dir, posre_file), self.restraints["posres"]
                )
            elif len(posre_strings.intersection(set(topology.binder_chains))) != 0:
                self.__overwrite_posres_file__(
                    Path(self.dir, posre_file), self.restraints["posres"]
                )
            else:
                self.__overwrite_posres_file__(
                    Path(self.dir, posre_file), self.restraints["posres_water"]
                )
        zip_top_path = Path(self.dir, f"restrained_{self.out_name}.zip")
        with ZipFile(zip_top_path, mode="w") as zf:
            for posre_file in filelist:
                zf.write(Path(self.dir, posre_file), arcname=posre_file)
        return ZipTopology.from_path_with_chains(
            zip_top_path,
            target_chains=topology.target_chains,
            binder_chains=topology.binder_chains,
        )

    @staticmethod
    def __overwrite_posres_file__(posre_file: Path, new_restraint: float) -> None:
        if new_restraint != 1000:
            set_posres(
                posre_file, Path(posre_file.parent, "temp_posres.itp"), new_restraint
            )
            shutil.move(Path(posre_file.parent, "temp_posres.itp"), posre_file)
        return


@define(frozen=True)
class BranchMDParams:
    gpu_id: int = field(converter=int, kw_only=True)
    omp_threads: int = field(converter=int, kw_only=True)
    mpi_threads: int = field(converter=int, kw_only=True)
    pinoffset: int = field(converter=int, kw_only=True)


def get_ngpus():
    try:
        p = sp.run(
            "nvidia-smi -L | wc -l",
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            shell=True,
            text=True,
        )
        ngpus = int(p.stdout)
        assert ngpus != 0
    except (ValueError, AssertionError, Exception):
        raise RuntimeError("No GPUs detected. Can't run locuaz.")
    return ngpus


def round_to_element_right(v: List[int], n: int) -> int:
    idx = np.searchsorted(v, n, side="right")
    try:
        return v[idx]
    except IndexError:
        return v[-1]


def numa_partition(
    all_threads, ngpus: int, nbranches: int, numa_regions: int
) -> Tuple[List[int], int, int]:
    numa_threads = all_threads[0 : len(all_threads) // numa_regions]
    numa_ngpus = ceil(ngpus / numa_regions)
    numa_nbranches = ceil(nbranches / numa_regions)

    return numa_threads, numa_ngpus, numa_nbranches


def auto_md_params(
    config_md: Dict[str, any],
    all_threads: List[int],
    nbranches: int,
    ngpus: Optional[int] = None,
) -> Tuple[List[int], int, int, int, int]:
    numa_regions = config_md["numa_regions"]
    if not ngpus:
        ngpus = get_ngpus()
    numa_threads, numa_ngpus, numa_nbranches = numa_partition(
        all_threads, ngpus, nbranches, numa_regions
    )
    # Define the number of threads for each run and try to maximize the number of parallel runs.
    max_threads_per_branch = 2
    for nbranches_parallel in range(numa_nbranches, 0, -1):
        try:
            # Get the max number of OMP threads, while also giving 1 extra thread for each parallel run (for MPI).
            max_threads_per_branch = (
                len(numa_threads) - nbranches_parallel
            ) // nbranches_parallel
        except ZeroDivisionError:
            raise RuntimeError(
                "Need at least 2 threads to use as OMP_THREADS. This shouldn't happen."
            )
        if max_threads_per_branch > 1:
            break

    threads_per_branch = round_to_element_right(OMP_MPI_THREADS, max_threads_per_branch)
    # Leave 2 threads between each branch, for the MPI process
    omp_threads = threads_per_branch - 2
    # Get the actual number of branches we can run parallely on each GPU
    nbranches_per_numa = len(numa_threads) // threads_per_branch
    # Get the pinoffsets for a NUMA section
    numa_pin_offsets = np.array(
        numa_threads[0::threads_per_branch][0:nbranches_per_numa]
    )
    # And now all of them:
    numa_step = len(all_threads) // numa_regions
    pin_offsets = list(
        chain.from_iterable(
            [list(numa_pin_offsets + i * numa_step) for i in range(numa_regions)]
        )
    )
    # Get the actual number of branches that could be run in parallel.
    # ``nbranches_parallel`` may overestimate the actual number of branches that will be run in parallel.
    # Eg: ``nbranches=7``, ``gpus = 4``, ``nbranches_per_numa = 2`` will give a ``nbranches_parallel=8``.
    # The first 3 GPUs will be running 2 branches at the same time, and the fourth one just one.
    nbranches_parallel = numa_regions * nbranches_per_numa
    #  Finally, repeat the pin offsets in case we have to run more branches that we can run simultaneously.
    nbr_of_rounds = ceil(nbranches / nbranches_parallel)
    pin_offsets = pin_offsets * nbr_of_rounds

    return pin_offsets, ngpus, omp_threads, 1, nbranches_parallel


def get_md_params(
    config_md: Dict[str, any], epoch: Epoch
) -> Tuple[Dict[str, BranchMDParams], int]:
    nbranches = len(epoch)
    epoch_md_params = {}
    if config_md["mps"]:
        all_threads = sorted(list(os.sched_getaffinity(0)))
        (
            pin_offsets,
            ngpus,
            omp_threads,
            mpi_threads,
            max_parallel_workers,
        ) = auto_md_params(config_md, all_threads, nbranches)
        for idx, (branch_name, branch) in enumerate(epoch.items()):
            epoch_md_params[branch_name] = BranchMDParams(
                gpu_id=idx % ngpus,
                omp_threads=omp_threads,
                mpi_threads=mpi_threads,
                pinoffset=pin_offsets[idx],
            )
    else:
        max_parallel_workers = config_md["ngpus"]
        for idx, (branch_name, branch) in enumerate(epoch.items()):
            epoch_md_params[branch_name] = BranchMDParams(
                gpu_id=idx % config_md["ngpus"],
                omp_threads=config_md["omp_procs"],
                mpi_threads=config_md["mpi_procs"],
                pinoffset=config_md["pinoffsets"][idx % config_md["ngpus"]],
            )

    return epoch_md_params, max_parallel_workers
