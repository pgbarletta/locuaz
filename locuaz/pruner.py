from typing import Tuple, Dict
from queue import PriorityQueue
import logging

from locuaz.projectutils import Epoch, Iteration

##### DEPRECATED ######

def top_pruner(better_iters: PriorityQueue, prune: int) -> Dict[str, Iteration]:
    """top_pruner select, at least, the first `prun` iterations above the threshold

    Args:
        better_iters (PriorityQueue[Tuple[int, Iteration]]): iterations above the threshold,
            sorted by approved SFs
        prune (int): minimum number of iterations returned

    Returns:
        Dict[str, Iteration]: map of top `iter_name:Iteration` for the next epoch.
    """

    top_iterations: Dict[str, Iteration] = {}
    for i in range(prune):
        iter = better_iters.get()[1]
        top_iterations[iter.iter_name] = iter
    return top_iterations


def adaptive_pruner(better_iters: PriorityQueue, prune: int) -> Dict[str, Iteration]:
    """adaptive_pruner select all iterations that have the maximum number of approved SFs

    Args:
        better_iters (PriorityQueue[Tuple[int, Iteration]]): iterations above the threshold,
            sorted by approved SFs
        prune (int): unused

    Returns:
        Dict[str, Iteration]: map of top `iter_name:Iteration` for the next epoch.
    """
    top_iterations: Dict[str, Iteration] = {}
    prev_count = 1
    while not better_iters.empty():
        # Remember, `count`, the priority, is negative.
        count, iter = better_iters.get()
        if count > prev_count:
            assert len(top_iterations) > 0, f"Logical error. This can't happen."
            break
        prev_count = count
        top_iterations[iter.iter_name] = iter
    return top_iterations
