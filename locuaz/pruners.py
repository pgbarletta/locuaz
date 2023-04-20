from typing import Dict, Callable, Set
import logging

from locuaz.projectutils import WorkProject
from locuaz.pruner import choose_top_iters, adaptive_pruner, top_pruner, threshold_pruner


pruners: Dict[str, Callable] = {
    "adaptive": adaptive_pruner,
    "top": top_pruner,
    "threshold": threshold_pruner,
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
    failed_pos: Set[int] = set()
    if better_iters_queue.empty():
        # All new iterations are worse than the old ones or they all failed during MD.
        log.info(
            f"Failed epoch after mutating resSeqs: {this_epoch.mutated_positions}. "
            f"Backing up epoch {this_epoch.id}."
        )
        failed_pos = this_epoch.mutated_positions
        this_epoch.backup()
        work_pjct.epochs[-1] = prev_epoch
    else:
        log.info(
            f"Successful epoch after mutating resSeqs: {this_epoch.mutated_positions} ."
        )
        pruner_func = pruners[work_pjct.config["protocol"]["pruner"]]
        prune = work_pjct.config["protocol"].get("prune")
        # Prune as many branches as requested, using the required pruner.
        this_epoch.top_iterations = pruner_func(better_iters_queue, prune)

    if work_pjct.has_failed_memory:
        work_pjct.failed_mutated_positions.appendleft(failed_pos)

    top_itrs_str = " ; ".join(
        [
            f"{iter.epoch_id}-{iter.iter_name}"
            for iter in work_pjct.epochs[-1].top_iterations.values()
        ]
    )
    log.info(f"Top iterations: {top_itrs_str}")

    return
