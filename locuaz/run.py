import logging
from pathlib import Path
from multiprocessing import Queue, SimpleQueue
import queue
from concurrent.futures import ProcessPoolExecutor
import concurrent.futures as cf

from projectutils import Epoch, Iteration, WorkProject
from runutils import MDrun

def run_epoch(work_pjct: WorkProject) -> None:
    run_min_nvt_epoch(work_pjct)
    run_npt_epoch(work_pjct)


def worker(task_queue: Queue, results_queue: SimpleQueue):
    while True:
        try:
            task = task_queue.get(block=False)
        except queue.Empty:
            print("Queue is empty! My work here is done. Exiting.")
            return
        engine, input_complex = task["engine"], task["complex"]

        complex = engine(input_complex)
        results_queue.put(complex)
        task_queue.task_done()

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
        for idx, (iter_name, iter) in enumerate(epoch.items()):
            
            gpu_id[iter_name] = idx % ngpus
            pinoffset[iter_name] = pinoffsets[idx % ngpus]
            
            log.info(f"Queuing MIN of the iteration {iter_name} to GPU {gpu_id[iter_name]} "
            f"and pinoffset: {pinoffset[iter_name]}.")

            min = MDrun.min(iter.dir_handle, work_pjct=work_pjct, gpu_id = gpu_id[iter_name],
                pinoffset=pinoffset[iter_name], out_name="min_" + work_pjct.name)
            futuros_min.append(ex.submit(min, iter.complex))
    
        for futu_min in cf.as_completed(futuros_min):
            if futu_min.exception():
                log.error(f"Exception while running MIN: {futu_min.exception()}")
                
            min_complex = futu_min.result()
            iter_name = '-'.join(min_complex.dir.dir_path.name.split('-')[1:])
            iter = epoch[iter_name]
            
            log.info(f"Queuing NVT of the iteration {iter_name} to GPU {gpu_id[iter_name]} "
            f"and pinoffset: {pinoffset[iter_name]}.")

            nvt = MDrun.nvt(
                iter.dir_handle, work_pjct=work_pjct, gpu_id = gpu_id[iter_name],
                pinoffset=pinoffset[iter_name], out_name="nvt_" + work_pjct.name)
            futuros_nvt.append(ex.submit(nvt, min_complex))

        for futu_nvt in cf.as_completed(futuros_nvt):
            if futu_nvt.exception():
                log.error(f"Exception while running NVT:  {futu_nvt.exception()} ")
                
            nvt_complex = futu_nvt.result()
            iter_name = '-'.join(nvt_complex.dir.dir_path.name.split('-')[1:])
            epoch[iter_name].complex = nvt_complex


def run_npt_epoch(work_pjct: WorkProject) -> None:
    log = logging.getLogger(work_pjct.name)
    epoch = work_pjct.epochs[-1]
    ngpus = work_pjct.config["md"]["ngpus"]
    prefix = work_pjct.config["main"]["prefix"]
    pinoffsets = work_pjct.config["md"]["pinoffsets"]

    with ProcessPoolExecutor(max_workers=ngpus) as ex:
        gpu_id = {}
        pinoffset = {}
        futuros_npt = []
        for idx, (iter_name, iter) in enumerate(epoch.items()):
            if iter.npt_started:
                continue
            gpu_nbr = idx % ngpus
            gpu_id[iter_name] = idx % ngpus
            pinoffset[iter_name] = pinoffsets[idx % ngpus]
            log.info(f"Queuing NPT of the iteration {iter_name} to GPU {gpu_nbr}. "
            f"and pinoffset: {pinoffset[iter_name]}.")
            
            npt = MDrun.npt(
                iter.dir_handle, work_pjct=work_pjct, gpu_id = gpu_nbr,
                pinoffset=pinoffset[iter_name], out_name=prefix + work_pjct.name)
            futuros_npt.append(ex.submit(npt, iter.complex))

        for futu_npt in cf.as_completed(futuros_npt):
            if futu_npt.exception():
                log.error(f"Exception while running NPT:  {futu_npt.exception()} ")
                
            npt_complex = futu_npt.result()
            iter_name = '-'.join(npt_complex.dir.dir_path.name.split('-')[1:])
            iter = epoch[iter_name]
            iter.complex = npt_complex


        