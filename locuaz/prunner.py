from ast import Call
from typing import Tuple, Dict, Callable
from queue import PriorityQueue

from projectutils import Epoch, Iteration, WorkProject


def beats_old_iter(
    old_iter: Iteration, new_iter: Iteration, threshold: int
) -> Tuple[bool, int]:

    # Allow the user to change scoring functions mid-run and only use
    # the subset present in both iterations.
    scoring_functions = set(old_iter.scores.keys()) & set(new_iter.scores.keys())
    count = sum(
        [
            old_iter.mean_scores[SF] > new_iter.mean_scores[SF]
            for SF in scoring_functions
        ]
    )

    return count > threshold, count


def choose_top_iters(
    prev_epoch: Epoch, this_epoch: Epoch, threshold: int
) -> PriorityQueue:
    better_iters: PriorityQueue = PriorityQueue()

    for iter in this_epoch.values():
        better_overall: bool = True
        rank_overall: int = 0
        for prev_iter in prev_epoch.top_iterations.values():
            better, rank = beats_old_iter(prev_iter, iter, threshold)
            better_overall &= better
            rank_overall += rank
        if better_overall:
            better_iters.put((-rank_overall, iter))

    return better_iters


def top_prunner(better_iters: PriorityQueue, prune: int) -> Dict[str, Iteration]:
    """top_prunner select, at least, the first `prun` iterations above the threshold

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


def adaptive_prunner(better_iters: PriorityQueue, prune: int) -> Dict[str, Iteration]:
    """adaptive_prunner select all iterations that have the maximum number of approved SFs

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
