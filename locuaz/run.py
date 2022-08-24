from asyncio import as_completed
import logging

# from queue import Queue, SimpleQueue, Empty
from multiprocessing import Queue, SimpleQueue
import queue
from concurrent.futures import ProcessPoolExecutor
import concurrent.futures as cf

from projectutils import Epoch, Iteration, WorkProject
from runutils import MDrun


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


# def minimize_iterations(work_pjct: WorkProject):
#     epoch = work_pjct.epochs[-1]
#     ngpus = work_pjct.config["md"]["ngpus"]

#     work_to_do: Queue = Queue()
#     for iter_name, iter in epoch.items():
#         logging.info(f"Queing minimization of the iteration {iter_name}.")

#         min = MDrun.min(
#             iter.dir_handle, work_pjct=work_pjct, out_name="min_" + work_pjct.name
#         )
#         task = {"engine": min, "complex": iter.complex}
#         work_to_do.put(task)

#     work_done: SimpleQueue = SimpleQueue()
#     with ProcessPoolExecutor(max_workers=ngpus) as ex:

#         futuros = [ex.submit(worker, work_to_do, work_done) for _ in range(ngpus)]
#         for futu in cf.as_completed(futuros):
#             print(f"{futu}")

def minimize_iterations(work_pjct: WorkProject):
    epoch = work_pjct.epochs[-1]
    ngpus = work_pjct.config["md"]["ngpus"]

    with ProcessPoolExecutor(max_workers=ngpus) as ex:
        futuros = []
        for iter_name, iter in epoch.items():
            logging.info(f"Queing minimization of the iteration {iter_name}.")

            min = MDrun.min(
                iter.dir_handle, work_pjct=work_pjct, out_name="min_" + work_pjct.name
            )
            futuros.append(ex.submit(min, iter.complex))
    
        for futu in cf.as_completed(futuros):
            if futu.exception():
                logging.error(f"")
            print(f"{futu}")

        