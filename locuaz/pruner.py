from typing import Tuple, Dict
from queue import PriorityQueue
import logging

from locuaz.projectutils import Epoch, Branch

##### DEPRECATED ######

def top_pruner(better_branches: PriorityQueue, prune: int) -> Dict[str, Branch]:
    """top_pruner select, at least, the first `prun` branches above the threshold

    Args:
        better_branches (PriorityQueue[Tuple[int, Branch]]): branches above the threshold,
            sorted by approved SFs
        prune (int): minimum number of branches returned

    Returns:
        Dict[str, Branch]: map of top `branch_name:Branch` for the next epoch.
    """

    top_branches: Dict[str, Branch] = {}
    for i in range(prune):
        iter = better_branches.get()[1]
        top_branches[iter.branch_name] = iter
    return top_branches


def adaptive_pruner(better_branches: PriorityQueue, prune: int) -> Dict[str, Branch]:
    """adaptive_pruner select all branches that have the maximum number of approved SFs

    Args:
        better_branches (PriorityQueue[Tuple[int, Branch]]): branches above the threshold,
            sorted by approved SFs
        prune (int): unused

    Returns:
        Dict[str, Branch]: map of top `branch_name:Branch` for the next epoch.
    """
    top_branches: Dict[str, Branch] = {}
    prev_count = 1
    while not better_branches.empty():
        # Remember, `count`, the priority, is negative.
        count, iter = better_branches.get()
        if count > prev_count:
            assert len(top_branches) > 0, f"Logical error. This can't happen."
            break
        prev_count = count
        top_branches[iter.branch_name] = iter
    return top_branches
