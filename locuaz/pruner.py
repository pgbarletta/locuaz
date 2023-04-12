from typing import Tuple, Dict
from queue import PriorityQueue
import logging

from projectutils import Epoch, Iteration


def beats_old_iter(
    old_iter: Iteration, new_iter: Iteration, threshold: int, log: logging.Logger
) -> Tuple[bool, int]:

    # Allow the user to change scoring functions mid-run and only use
    # the subset present in both iterations.
    old_SFs = set(old_iter.scores.keys())
    new_SFs = set(new_iter.scores.keys())
    scoring_functions = old_SFs & new_SFs
    if len(scoring_functions) == 0:
        raise RuntimeError(f"No common scoring functions between the ones from the old iteration ({old_SFs}) "
                           f"and those from the new one ({new_SFs}). Cannot prune.")
    log.info(f"Scoring functions: {scoring_functions}")
    count = sum(
        [
            old_iter.mean_scores[SF] >= new_iter.mean_scores[SF]
            for SF in scoring_functions
        ]
    )
    log.info(
        f"{new_iter.epoch_id}-{new_iter.iter_name} vs. {old_iter.epoch_id}-{old_iter.iter_name} "
        f"improves on {count} of {len(scoring_functions)} scoring functions."
    )

    return count >= threshold, count


def choose_top_iters(
    prev_epoch: Epoch, this_epoch: Epoch, threshold: int, log: logging.Logger
) -> PriorityQueue:
    better_iters: PriorityQueue = PriorityQueue()

    for iter in this_epoch.values():
        better_overall: bool = True
        rank_overall: int = 0
        for prev_iter in prev_epoch.top_iterations.values():
            better, rank = beats_old_iter(prev_iter, iter, threshold, log)
            better_overall &= better
            rank_overall += rank
        if better_overall:
            better_iters.put((-rank_overall, iter))

    return better_iters


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


def threshold_pruner(better_iters: PriorityQueue, prune: int) -> Dict[str, Iteration]:
    """threshold_pruner select all iterations that satisfy the threshold of approved SFs

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
        _, iter = better_iters.get()
        top_iterations[iter.iter_name] = iter
    return top_iterations
