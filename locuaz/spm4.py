from collections.abc import Iterable, ItemsView
from typing import List, Tuple, Dict, Set, Iterator
from random import choice, sample
from collections import defaultdict
from logging import Logger
from deprecated.sphinx import deprecated

from locuaz.projectutils import Branch, Epoch
from locuaz.mutation import Mutation
from locuaz.abstractmutationgenerator import AbstractMutationGenerator

@deprecated(version="0.7.0",
            reason="Mutation Generators are replaced by Mutation Creators.")
class SPM4(AbstractMutationGenerator):
    """
    SPM4 generates mutations by splitting all amino acids (except CYS) in the following categories:
    negative, positive, hydrophobic and ring-containing.
    It will pick 1 amino acid from each category until it generates one for each branch and randomly pick a
    position for all of them, hence SPM (Single Point Mutation).
    """
    NEG_AAS: Tuple = ("D", "E", "S", "T")
    POS_AAS: Tuple = ("R", "N", "Q", "H", "K")
    # Cystein is not included.
    PHO_AAS: Tuple = ("A", "G", "I", "M", "L", "V")
    RIN_AAS: Tuple = ("P", "F", "W", "Y")
    CAT_AAS: Tuple[Tuple, Tuple, Tuple, Tuple] = (
        NEG_AAS,
        POS_AAS,
        PHO_AAS,
        RIN_AAS,
    )
    N_CAT = len(CAT_AAS)
    remaining_categories: Set[int]
    excluded_aas: Set[str]
    excluded_pos: Set[int]
    idx_chain: int
    idx_residue: int

    mutations: Dict[str, List[Mutation]]

    def __init__(
            self,
            epoch: Epoch,
            branches: int,
            *,
            excluded_aas: Set[str],
            excluded_pos: Set[int],
            probe_radius: float = 1.4,
            use_tleap: bool = False,
            logger: Logger = None
    ) -> None:
        self.remaining_categories = set(range(0, self.N_CAT))
        self.mutations = defaultdict(list)
        self.excluded_aas = excluded_aas
        self.excluded_pos = excluded_pos

        mut_idx_chain, mut_chainID, mut_idx_residue, mut_resSeq = self.__generate_position__(epoch, use_tleap, logger)
        self.__fill_mutations__(epoch, branches, mut_idx_chain, mut_chainID, mut_idx_residue, mut_resSeq)

    def __fill_mutations__(self, epoch: Epoch, branches: int, mut_idx_chain: int, mut_chainID: str,
                           mut_idx_residue: int, mut_resSeq: int):
        """
        generates up to `branches` different mutations
        Parameters
        ----------
        epoch :
        branches :
        mut_idx_chain :
        mut_chainID :
        mut_idx_residue :
        mut_resSeq :

        Returns
        -------

        """
        new_branches = branches
        remaining_branches = set(epoch.top_branches.keys())
        while new_branches != 0:
            branch = epoch.top_branches[remaining_branches.pop()]
            old_aa, new_aa = self.__get_random_aa__(branch, mut_idx_chain, mut_idx_residue)
            # Build the mutation object for this branch.
            mutation = Mutation(
                chainID=mut_chainID,
                resSeq=mut_resSeq,
                old_aa=old_aa,
                new_aa=new_aa,
                chainID_idx=mut_idx_chain,
                resSeq_idx=mut_idx_residue,
            )

            self.mutations[branch.branch_name].append(mutation)
            if len(remaining_branches) == 0:
                # If all branches have been mutated at least once, and we still
                # have branches to generate, restart `remaining_branches`.
                remaining_branches = set(epoch.top_branches.keys())
            new_branches -= 1

    def __generate_position__(self, epoch: Epoch, use_tleap: bool = False,
                              logger: Logger = None) -> Tuple[int, str, int, int]:
        # Get a branch to read the chainIDs and the resSeqs.
        try:
            any_branch = next(iter(epoch.top_branches.values()))
        except Exception:
            raise RuntimeError(f"No available branches on Epoch {epoch.id}. "
                               "It's likely that all of them failed during MD.")

        # Now, filter the mutating resSeqs.
        candidates_resSeq = [resSeq for cdr in any_branch.resSeqs for resSeq in cdr if
                             resSeq not in self.excluded_pos]
        if len(candidates_resSeq) == 0:
            raise RuntimeError(f"Cannot mutate. No CDR residue that isn't excluded.")

        logger.info(f"Generating mutations with: {self}.\n"
                    f"'mutating_resSeq': {any_branch.resSeqs}.\n"
                    f"excluded resSeq: {self.excluded_pos}.\n"
                    f"'mutating_resSeq' that may be mutated: {candidates_resSeq}.")

        # Choose the position to mutate. This will be the same for all branches.
        mut_resSeq = choice(candidates_resSeq)
        # Now, get the remaining details associated to `mut_resSeq`, including the chain from where it came from
        for j, resSeqs in enumerate(any_branch.resSeqs):
            if mut_resSeq in resSeqs:
                mut_idx_chain = j
                mut_chainID = any_branch.chainIDs[mut_idx_chain]
                mut_idx_residue = resSeqs.index(mut_resSeq)
                return mut_idx_chain, mut_chainID, mut_idx_residue, mut_resSeq

    def __pop_random_aa__(self) -> str:
        cat_idx = choice(tuple(self.remaining_categories))
        self.remaining_categories.difference_update({cat_idx})

        shuffled_aas = sample(self.CAT_AAS[cat_idx], len(self.CAT_AAS[cat_idx]))
        for new_aa in shuffled_aas:
            if new_aa not in self.excluded_aas:
                self.excluded_aas.add(new_aa)
                return new_aa
        raise RuntimeError("Can't generate novel AA. This shouldn't happen.")

    def __get_random_aa__(self, branch: Branch, idx_chain: int, idx_residue: int) -> Tuple[str, str]:
        old_aa = branch.resnames[idx_chain][idx_residue]
        self.excluded_aas.add(old_aa)
        new_aa = self.__pop_random_aa__()

        if len(self.remaining_categories) == 0:
            # All categories have already been chosen from `N` times.
            # Allow all of them again for the `N+1` branch, except those
            # that are exhausted already.
            self.remaining_categories = set(range(self.N_CAT))
            for i in range(self.N_CAT):
                if set(self.CAT_AAS[i]).issubset(self.excluded_aas):
                    # All AAs from this category have already been chosen
                    self.remaining_categories.difference_update({i})
        return old_aa, new_aa

    def __getitem__(self, key: str) -> Iterable[Mutation]:
        return self.mutations[key]

    def __iter__(self) -> Iterator:
        return self.mutations.__iter__()

    def __contains__(self, value: Branch) -> bool:
        return self.mutations.__contains__(value)

    def __len__(self) -> int:
        return len(self.mutations)

    def __items__(self) -> ItemsView:
        return self.mutations.items()

    def __str__(self) -> str:
        return "Single Point Mutation from 4 amino acid groups"
