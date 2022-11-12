from typing import Dict, Callable
import logging

from projectutils import WorkProject
from prunner import choose_top_iters, adaptive_prunner, top_prunner


prunners: Dict[str, Callable] = {
    "adaptive": adaptive_prunner,
    "top": top_prunner,
}


def prune(work_pjct: WorkProject) -> None:
    this_epoch = work_pjct.epochs[-1]
    try:
        prev_epoch = work_pjct.epochs[-2]
    except IndexError as e:
        # Initial epoch.
        this_epoch.top_iterations = this_epoch.iterations
        return

    log = logging.getLogger(work_pjct.name)
    threshold = work_pjct.config["scoring"]["consensus_threshold"]

    better_iters_queue = choose_top_iters(prev_epoch, this_epoch, threshold, log)
    if better_iters_queue.empty():
        # All new iterations are worse than the old ones or they all failed during MD.
        log.info(f"Failed epoch. Backing up epoch {this_epoch.id}.")
        prev_epoch.backup()
        this_epoch = prev_epoch
    else:
        prunner_func = prunners[work_pjct.config["protocol"]["prunner"]]
        prune = work_pjct.config["protocol"].get("prune")
        # Use the required prunner
        this_epoch.top_iterations = prunner_func(better_iters_queue, prune)

    top_itrs_str = " ; ".join(
        [
            f"{iter.epoch_id}-{iter.iter_name}"
            for iter in this_epoch.top_iterations.values()
        ]
    )
    log.info(f"Top iterations: {top_itrs_str}")

    return
