from typing import Set, Dict
import logging

from locuaz.projectutils import WorkProject, Branch
from locuaz.pruners import pruners


def prune(work_pjct: WorkProject) -> None:
    this_epoch = work_pjct.epochs[-1]
    try:
        prev_epoch = work_pjct.epochs[-2]
    except IndexError:
        # Initial epoch.
        this_epoch.top_branches = this_epoch.branches
        return

    log = logging.getLogger(work_pjct.name)
    # Use chosen pruner to prune.
    pruner_func = pruners[work_pjct.config["pruning"]["pruner"]](prev_epoch, this_epoch, log)
    better_branches_queue = pruner_func.prune(work_pjct.config)

    failed_pos: Set[int] = set()
    if better_branches_queue.empty():
        # All new branches are worse than the old ones, or they all failed during MD.
        log.info(f"Unsuccessful epoch after mutating resSeqs: {this_epoch.mutated_positions}. "
                 f"Backing up epoch {this_epoch.id}.")
        failed_pos = this_epoch.mutated_positions
        this_epoch.backup()
        work_pjct.epochs[-1] = prev_epoch
    else:
        log.info(f"Successful epoch after mutating resSeqs: {this_epoch.mutated_positions} .")
        this_epoch.top_branches: Dict[str, Branch] = {}
        while not better_branches_queue.empty():
            _, branch = better_branches_queue.get()
            this_epoch.top_branches[branch.branch_name] = branch

    if work_pjct.has_failed_memory:
        # If epoch was successful, we're appending an empty set, just to move the queue along.
        work_pjct.failed_mutated_positions.appendleft(failed_pos)

    top_itrs_str = " ; ".join([f"{branch.epoch_id}-{branch.branch_name}"
                               for branch in work_pjct.epochs[-1].top_branches.values()])
    log.info(f"Top branches: {top_itrs_str}")

    return
