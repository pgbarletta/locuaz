from typing import Set, Dict
import logging

from locuaz.projectutils import WorkProject, Iteration
from locuaz.pruners import pruners


def prune(work_pjct: WorkProject) -> None:
    this_epoch = work_pjct.epochs[-1]
    try:
        prev_epoch = work_pjct.epochs[-2]
    except IndexError:
        # Initial epoch.
        this_epoch.top_iterations = this_epoch.iterations
        return

    log = logging.getLogger(work_pjct.name)
    # Use chosen pruner to prune.
    pruner_func = pruners[work_pjct.config["pruning"]["pruner"]](prev_epoch, this_epoch, log)
    better_iters_queue = pruner_func.prune(work_pjct.config)

    failed_pos: Set[int] = set()
    if better_iters_queue.empty():
        # All new iterations are worse than the old ones, or they all failed during MD.
        log.info(f"Unsuccessful epoch after mutating resSeqs: {this_epoch.mutated_positions}. "
                 f"Backing up epoch {this_epoch.id}.")
        failed_pos = this_epoch.mutated_positions
        this_epoch.backup()
        work_pjct.epochs[-1] = prev_epoch
    else:
        log.info(f"Successful epoch after mutating resSeqs: {this_epoch.mutated_positions} .")
        this_epoch.top_iterations: Dict[str, Iteration] = {}
        while not better_iters_queue.empty():
            _, iteration = better_iters_queue.get()
            this_epoch.top_iterations[iteration.iter_name] = iteration

    if work_pjct.has_failed_memory:
        # If epoch was successful, we're appending an empty set, just to move the queue along.
        work_pjct.failed_mutated_positions.appendleft(failed_pos)

    top_itrs_str = " ; ".join([f"{iteration.epoch_id}-{iteration.iter_name}"
                               for iteration in work_pjct.epochs[-1].top_iterations.values()])
    log.info(f"Top iterations: {top_itrs_str}")

    return