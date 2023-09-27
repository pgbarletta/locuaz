import logging
from pathlib import Path
from typing import Dict, Any, Optional
from concurrent.futures import ProcessPoolExecutor
import concurrent.futures as cf

from locuaz.projectutils import WorkProject, Epoch
from locuaz.runutils import MDrun, BranchMDParams, get_md_params


def run_epoch(work_pjct: WorkProject) -> None:
    log = logging.getLogger(f"{work_pjct.name}")
    epoch = work_pjct.epochs[-1]
    md_params, max_parallel_workers = get_md_params(work_pjct.config["md"], epoch)
    if not epoch.nvt_done:
        run_min_nvt_epoch(epoch, work_pjct.config["md"], md_params, max_parallel_workers, work_pjct.name, log)
        run_npt_epoch(epoch, work_pjct.config["md"], md_params, max_parallel_workers, work_pjct.name, log)
    elif not epoch.npt_done:
        log.info(f"Skipping minimization and NVT run of epoch {epoch.id}")
        run_npt_epoch(epoch, work_pjct.config["md"], md_params, max_parallel_workers, work_pjct.name, log)
    else:
        log.info(f"Skipping NPT run of epoch {epoch.id}")


def run_min_nvt_epoch(epoch: Epoch, md_config: Dict[str, Any], md_params: Dict[str, BranchMDParams],
                      max_parallel_workers: int, name: str, log: Optional[logging.Logger] = None) -> None:
    with ProcessPoolExecutor(max_workers=max_parallel_workers) as ex:
        futuros_min = []
        futuros_nvt = []
        for idx, (branch_name, branch) in enumerate(epoch.items()):
            gpu_id = md_params[branch_name].gpu_id
            omp_threads = md_params[branch_name].omp_threads
            mpi_threads = md_params[branch_name].mpi_threads
            pinoffset = md_params[branch_name].pinoffset
            # Ideally the MDrun class would do the logging, but python's logging is not thread safe.
            if log:
                log.info(f"Queuing MIN of the branch {branch_name}. "
                         f"{gpu_id=}, {omp_threads=}, {mpi_threads=}, {pinoffset=}.")

            min = MDrun.min(branch.dir_handle, gmx_mdrun=md_config["gmx_mdrun"],
                            min_mdp=md_config["mdp_paths"]["min_mdp"], gpu_id=gpu_id, omp_threads=omp_threads,
                            mpi_threads=mpi_threads, pinoffset=pinoffset, out_name=f"min_{name}",
                            maxwarn=md_config["maxwarn"])
            futuros_min.append(ex.submit(min, branch.complex))

        for futu_min in cf.as_completed(futuros_min):
            if futu_min.exception():
                log.error(f"Exception while running MIN: {futu_min.exception()}")

            _, min_complex = futu_min.result()
            branch_name = '-'.join(min_complex.dir.dir_path.name.split('-')[1:])
            gpu_id = md_params[branch_name].gpu_id
            omp_threads = md_params[branch_name].omp_threads
            mpi_threads = md_params[branch_name].mpi_threads
            pinoffset = md_params[branch_name].pinoffset
            branch = epoch[branch_name]

            if log:
                log.info(
                    f"Queuing NVT of the branch {branch_name}. "
                    f"{gpu_id=}, {omp_threads=}, {mpi_threads=}, {pinoffset=}.")

            nvt = MDrun.nvt(branch.dir_handle, gmx_mdrun=md_config["gmx_mdrun"],
                            nvt_mdp=md_config["mdp_paths"]["nvt_mdp"], gpu_id=gpu_id, omp_threads=omp_threads,
                            mpi_threads=mpi_threads, pinoffset=pinoffset, out_name=f"nvt_{name}",
                            maxwarn=md_config["maxwarn"])
            futuros_nvt.append(ex.submit(nvt, min_complex))

            for futu_nvt in cf.as_completed(futuros_nvt):
                if futu_nvt.exception():
                    log.error(f"Exception while running NVT:  {futu_nvt.exception()} ")

            _, nvt_complex = futu_nvt.result()
            branch_name = '-'.join(nvt_complex.dir.dir_path.name.split('-')[1:])
            epoch[branch_name].complex = nvt_complex
            epoch.nvt_done = True


def run_npt_epoch(epoch: Epoch, md_config: Dict[str, Any], md_params: Dict[str, BranchMDParams],
                  max_parallel_workers: int, name: str, log: Optional[logging.Logger] = None) -> None:
    box_type = md_config["box_type"]
    with ProcessPoolExecutor(max_workers=max_parallel_workers) as ex:
        futuros_npt = {}
        for idx, (branch_name, branch) in enumerate(epoch.items()):
            gpu_id = md_params[branch_name].gpu_id
            omp_threads = md_params[branch_name].omp_threads
            mpi_threads = md_params[branch_name].mpi_threads
            pinoffset = md_params[branch_name].pinoffset
            # Ideally the MDrun class would do the logging, but python's logging is not thread safe.
            if log:
                log.info(f"Queuing NPT of the branch {branch_name}. "
                         f"{gpu_id=}, {omp_threads=}, {mpi_threads=}, {pinoffset=}.")

            npt = MDrun.npt(branch.dir_handle, gmx_mdrun=md_config["gmx_mdrun"],
                            npt_mdp=md_config["mdp_paths"]["npt_mdp"], gpu_id=gpu_id, omp_threads=omp_threads,
                            mpi_threads=mpi_threads, pinoffset=pinoffset, out_name=f"npt_{name}",
                            maxwarn=md_config["maxwarn"], restraints=md_config.get("npt_restraints"))
            futu_npt = ex.submit(npt, branch.complex)
            futuros_npt[futu_npt] = branch_name

        for futu_npt in cf.as_completed(futuros_npt):
            branch_name = futuros_npt[futu_npt]
            if futu_npt.exception():
                log.error(f"Exception while running NPT from branch: {branch_name}\n{futu_npt.exception()}")
                del epoch[branch_name]
                continue
            all_atoms_in_box, npt_complex = futu_npt.result()

            branch = epoch[branch_name]
            branch.complex = npt_complex
            if not all_atoms_in_box and box_type == "triclinic":
                log.error(f"{epoch.id}-{branch_name} has atoms outside the box. This run will not be scored.")
                branch.outside_box = True

    epoch.npt_done = True
