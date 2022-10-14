from typing import Dict, Callable

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

    threshold = work_pjct.config["scoring"]["consensus_threshold"]

    better_iters_queue = choose_top_iters(prev_epoch, this_epoch, threshold)
    if better_iters_queue.empty():
        # All new iterations are worse than the old ones or they all failed during MD.
        top_iterations = prev_epoch.top_iterations
    else:
        prunner_func = prunners[work_pjct.config["protocol"]["prunner"]]
        prune = work_pjct.config["protocol"]["prune"]
        # Use the required prunner
        top_iterations = prunner_func(better_iters_queue, prune)

    this_epoch.top_iterations = top_iterations
    return
