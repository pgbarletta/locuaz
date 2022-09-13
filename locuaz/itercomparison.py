import logging
from pathlib import Path
from typing import List, Sequence, Set, Dict, Tuple, Optional
from collections.abc import Sequence
from queue import PriorityQueue

from projectutils import Epoch, Iteration, WorkProject
from fileutils import DirHandle, FileHandle, copy_to, catenate, update_header


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


def choose_top_iters(work_pjct: WorkProject) -> None:

    this_epoch = work_pjct.epochs[-1]
    try:
        prev_epoch = work_pjct.epochs[-2]
    except IndexError as e:
        # Initial epoch.
        this_epoch.top_iterations = tuple(this_epoch.iterations.keys())
        return

    threshold = work_pjct.config["scoring"]["consensus_threshold"]
    prune = work_pjct.config["protocol"]["prune"]

    better_iters: PriorityQueue = PriorityQueue()
    for prev_iter_name in prev_epoch.top_iterations:
        prev_iter = prev_epoch[prev_iter_name]

        for iter_name, iter in this_epoch.items():
            better, sf_count = beats_old_iter(prev_iter, iter, threshold)
            if better:
                better_iters.put((sf_count, iter_name))

    if better_iters.empty():
        # All new iterations are worse than the old ones,
        # or they all failed during MD.
        this_epoch.top_iterations = prev_epoch.top_iterations
        this_epoch.iterations = prev_epoch.iterations
    else:
        this_epoch.top_iterations = tuple([better_iters.get()[1] for i in range(prune)])
