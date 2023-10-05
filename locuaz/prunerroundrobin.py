from queue import PriorityQueue
from typing import Dict
import logging

from locuaz.abstractpruner import AbstractPruner
from locuaz.projectutils import Epoch, Branch
from locuaz.primitives import PriorityDeque


class PrunerRoundRobin(AbstractPruner):
    """PrunerRoundRobin selects the top N branches after comparing the average scores
     of 'all_branches' and 'top branches' from the previous epoch in a round-robin
     approach (all against all). The selection of top N branches may include a
     mix of branches from the current epoch and the previous one.
    """

    def __init__(self, prev_epoch: Epoch, this_epoch: Epoch, log: logging.Logger):
        self.prev_epoch = prev_epoch
        self.this_epoch = this_epoch
        self.log = log

    def prune(self, config: Dict) -> PriorityQueue:
        """
        The only public method from the Pruner.

        Parameters
        ----------
        config : Dict
            User input config file

        Returns
        -------
        passing_branches: PriorityQueue
            ordered queue with the top branches from the new and the previous
            epochs. It may not be empty.
        """
        return self.__get_passing_branches__(config["pruning"]["roundrobin_threshold"])

    def __get_passing_branches__(self, roundrobin_threshold: int) -> PriorityQueue:
        """

        Parameters
        ----------
        roundrobin_threshold : int
            number of branches that will be selected as top branches.

        Returns
        -------
        passing_branches: PriorityQueue
            ordered queue with the branches from the new epoch that are better than all the
            branches from the old epoch. It may be empty.
        """
        better_branches: PriorityDeque = PriorityDeque(roundrobin_threshold)
        all_branches: Dict[str, Branch] = self.prev_epoch.top_branches | self.this_epoch
        for left_name, left_branch in all_branches.items():
            rank_overall: int = 0
            for right_name, right_branch in all_branches.items():
                if left_name == right_name:
                    continue
                rank = self.__beats_old_branch__(right_branch, left_branch)
                rank_overall += rank
            better_branches.put((-rank_overall, left_branch))

        # Cast it into PriorityQueue, since the code expects that
        passing_branches = PriorityQueue()
        [passing_branches.put(branch) for branch in better_branches]

        return passing_branches

    def __beats_old_branch__(self, old_branch: Branch, new_branch: Branch) -> int:
        """
        Parameters
        ----------
        old_branch: Branch
        new_branch : Branch

        Returns
        -------
        count: int
            number of improved scorers.
        """
        # Allow the user to change scorers mid-run and only use
        # the subset present in both branches.
        old_SFs = set(old_branch.scores.keys())
        new_SFs = set(new_branch.scores.keys())
        scorers = old_SFs & new_SFs
        if len(scorers) == 0:
            raise RuntimeError(f"No common scorers between the ones from the old branch ({old_SFs}) "
                               f"and those from the new one ({new_SFs}). Cannot prune.")
        self.log.info(f"Scorers {scorers}")
        count = sum([old_branch.mean_scores[SF] >= new_branch.mean_scores[SF]
                     for SF in scorers])
        self.log.info(
            f"{new_branch.epoch_id}-{new_branch.branch_name} vs. {old_branch.epoch_id}-{old_branch.branch_name} "
            f"improves on {count} of {len(scorers)} scorers.")

        return count
