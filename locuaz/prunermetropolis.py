from queue import PriorityQueue
from typing import Dict, Tuple
import logging
import numpy as np

from locuaz.abstractpruner import AbstractPruner
from locuaz.projectutils import Epoch, Branch


class PrunerMetropolis(AbstractPruner):
    """PrunerMetropolis uses only 1 scoring function (SF) and approves a branch when its score is better
    than the scores of all the previous ones, or in case it's not, if the metropolis acceptance criteria
    approves it.
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
            ordered queue with the branches from the new epoch that are better than all the
            branches from the old epoch. It may be empty.
        """
        return self.__get_passing_branches__(config["pruning"]["kT"])

    def __get_passing_branches__(self, kT: float = 0.593) -> PriorityQueue:
        """

        Parameters
        ----------
        kT : float
            boltzmann constant * temperature in Kcal/mol. Defaults to 0.593

        Returns
        -------
        passing_branches: PriorityQueue
            ordered queue with the branches from the new epoch that are better than all the
            branches from the old epoch. It may be empty.
        """
        better_branches: PriorityQueue = PriorityQueue()

        for branch in self.this_epoch.values():
            better_overall: bool = True
            rank_overall: int = 0
            for prev_branch in self.prev_epoch.top_branches.values():
                better, rank = self.__beats_old_branch__(prev_branch, branch, kT)
                better_overall &= better
                rank_overall += rank
            if better_overall:
                better_branches.put((-rank_overall, branch))

        return better_branches

    def __beats_old_branch__(self, old_branch: Branch, new_branch: Branch, kT: float) -> Tuple[bool, int]:
        """
        Parameters
        ----------
        old_branch: Branch
        new_branch : Branch
        kT : float
            boltzmann constant * temperature in Kcal/mol. Defaults to 0.593

        Returns
        -------
        beats and count: Tuple[bool, int]
            whether if it does beat the old iter and the number of improved SFs.
        """
        # Allow the user to change scorers mid-run and only use
        # the subset present in both branches.
        old_SFs = set(old_branch.scores.keys())
        new_SFs = set(new_branch.scores.keys())
        scorers = old_SFs & new_SFs
        if len(scorers) == 0:
            raise RuntimeError(f"No common scorers between the ones from the old branch ({old_SFs}) "
                               f"and those from the new one ({new_SFs}). Cannot prune.")
        self.log.info(f"Scoring function: {scorers}")

        sf = scorers.pop()
        old_energy = old_branch.mean_scores[sf]
        new_energy = new_branch.mean_scores[sf]
        if new_energy < old_energy:
            self.log.info(f"{new_branch.epoch_id}-{new_branch.branch_name} vs. {old_branch.epoch_id}-{old_branch.branch_name} "
                          f"improves on 1 scoring function.")
            return True, 1
        else:
            acceptance_probability = np.exp((old_energy - new_energy) / kT)
            if np.random.rand() < acceptance_probability:
                self.log.info(f"{new_branch.epoch_id}-{new_branch.branch_name} vs. {old_branch.epoch_id}-{old_branch.branch_name} "
                              f"is accepted.")
                return True, 0

        self.log.info(f"{new_branch.epoch_id}-{new_branch.branch_name} vs. {old_branch.epoch_id}-{old_branch.branch_name} "
                      f"is rejected.")
        return False, 0
