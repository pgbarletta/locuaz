import logging
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import concurrent.futures as cf

from locuaz.projectutils import WorkProject
from locuaz.runutils import MDrun

def run_epoch(work_pjct: WorkProject) -> None:
    log = logging.getLogger(f"{work_pjct.name}")
    
    if not work_pjct.epochs[-1].nvt_done:
        run_min_nvt_epoch(work_pjct)
    else:
        log.info(f"Skipping minimization and NVT run of epoch {work_pjct.epochs[-1].id}")
    
    if not work_pjct.epochs[-1].npt_done:
        run_npt_epoch(work_pjct)
    else:
        log.info(f"Skipping NPT run of epoch {work_pjct.epochs[-1].id}")

def run_min_nvt_epoch(work_pjct: WorkProject) -> None:
    log = logging.getLogger(f"{work_pjct.name}")
    epoch = work_pjct.epochs[-1]
    ngpus = work_pjct.config["md"]["ngpus"]
    pinoffsets = work_pjct.config["md"]["pinoffsets"]

    with ProcessPoolExecutor(max_workers=ngpus) as ex:
        futuros_min = []
        futuros_nvt = []
        gpu_id = {}
        pinoffset = {}
        for idx, (branch_name, iter) in enumerate(epoch.items()):
            
            gpu_id[branch_name] = idx % ngpus
            pinoffset[branch_name] = pinoffsets[idx % ngpus]
            
            log.info(f"Queuing MIN of the branch {branch_name} to GPU {gpu_id[branch_name]} "
            f"and pinoffset: {pinoffset[branch_name]}.")

            min = MDrun.min(iter.dir_handle, work_pjct=work_pjct, gpu_id = gpu_id[branch_name],
                pinoffset=pinoffset[branch_name], out_name="min_" + work_pjct.name)
            futuros_min.append(ex.submit(min, iter.complex))
    
        for futu_min in cf.as_completed(futuros_min):
            if futu_min.exception():
                log.error(f"Exception while running MIN: {futu_min.exception()}")
                
            _, min_complex = futu_min.result()
            branch_name = '-'.join(min_complex.dir.dir_path.name.split('-')[1:])
            iter = epoch[branch_name]
            
            log.info(f"Queuing NVT of the branch {branch_name} to GPU {gpu_id[branch_name]} "
            f"and pinoffset: {pinoffset[branch_name]}.")

            nvt = MDrun.nvt(
                Path(iter.dir_handle), work_pjct=work_pjct, gpu_id = gpu_id[branch_name],
                pinoffset=pinoffset[branch_name], out_name="nvt_" + work_pjct.name)
            futuros_nvt.append(ex.submit(nvt, min_complex))

        for futu_nvt in cf.as_completed(futuros_nvt):
            if futu_nvt.exception():
                log.error(f"Exception while running NVT:  {futu_nvt.exception()} ")
                
            _, nvt_complex = futu_nvt.result()
            branch_name = '-'.join(nvt_complex.dir.dir_path.name.split('-')[1:])
            epoch[branch_name].complex = nvt_complex
    epoch.nvt_done = True


def run_npt_epoch(work_pjct: WorkProject) -> None:
    log = logging.getLogger(work_pjct.name)
    epoch = work_pjct.epochs[-1]
    ngpus = work_pjct.config["md"]["ngpus"]
    prefix = work_pjct.config["main"]["prefix"]
    pinoffsets = work_pjct.config["md"]["pinoffsets"]
    box_type = work_pjct.config["md"]["box_type"]

    with ProcessPoolExecutor(max_workers=ngpus) as ex:
        gpu_id = {}
        pinoffset = {}
        futuros_npt = {}
        for idx, (branch_name, iter) in enumerate(epoch.items()):
            gpu_nbr = idx % ngpus
            gpu_id[branch_name] = idx % ngpus
            pinoffset[branch_name] = pinoffsets[idx % ngpus]
            log.info(f"Queuing NPT of the branch {branch_name} to GPU {gpu_nbr} "
            f"and pinoffset: {pinoffset[branch_name]}.")
            
            npt = MDrun.npt(
                iter.dir_handle, work_pjct=work_pjct, gpu_id = gpu_nbr,
                pinoffset=pinoffset[branch_name], out_name=prefix + work_pjct.name)
            
            futu_npt = ex.submit(npt, iter.complex)
            futuros_npt[futu_npt] = branch_name

        for futu_npt in cf.as_completed(futuros_npt):
            branch_name = futuros_npt[futu_npt]
            if futu_npt.exception():
                log.error(f"Exception while running NPT from branch: {branch_name}\n"
                    f"{futu_npt.exception()}")
                del epoch[branch_name]
                continue
            all_atoms_in_box, npt_complex = futu_npt.result()
            
            iter = epoch[branch_name]
            iter.complex = npt_complex
            if not all_atoms_in_box and box_type == "triclinic":
                log.error(f"{epoch.id}-{branch_name} has atoms outside the box. "
                "This run may not be apt to continue.")
                iter.outside_box = True

            
    epoch.npt_done = True


        
