from abc import ABC, abstractmethod
from email import iterators
from typing import List, Tuple, Dict, Set, Iterator, Sequence
from random import choice, sample
from collections import defaultdict
from collections.abc import ItemsView, Mapping

from projectutils import Iteration, Epoch
from mutator import Mutation


class AbstractMutationGenerator(ABC, Mapping):
    @abstractmethod
    def __init__(
        self,
        epoch: Epoch,
        max_branches: int,
        *,
        excluded_aas: Set[str],
        excluded_pos: Set[int],
    ) -> None:
        pass

    @abstractmethod
    def __getitem__(self, key: Iteration) -> Mutation:
        pass

    @abstractmethod
    def __iter__(self) -> Iterator:
        pass

    @abstractmethod
    def __contains__(self, value: Iteration) -> bool:
        pass


class SPM_RB(AbstractMutationGenerator):
    pass
    # def __init__(
    #     self,
    #     epoch: Epoch,
    #     max_branches: int,
    #     *,
    #     excluded_aas: Set[str],
    #     excluded_pos: Set[int],
    # ) -> None:
    #     pass

    # def __getitem__(self, key: Iteration) -> Mutation:
    #     pass

    # def __iter__(self) -> Iterator:
    #     pass

    # def __contains__(self, value: Iteration) -> bool:
    #     pass


class SPM_4(AbstractMutationGenerator):
    # Cystein is not included.
    NEG_AAS: Tuple = ("D", "E", "S", "T")
    POS_AAS: Tuple = ("R", "N", "Q", "H", "K")
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
        max_branches: int,
        *,
        excluded_aas: Set[str],
        excluded_pos: Set[int],
    ) -> None:
        self.remaining_categories = set(range(0, self.N_CAT))
        self.mutations = defaultdict(list)
        self.excluded_aas = excluded_aas
        self.excluded_pos = excluded_pos

        mut_chainID, mut_resSeq = self.__get_random_pos__(epoch)
        # Now, generate up to `max_branches` mutations
        remaining_branches = max_branches
        remaining_iterations = set(epoch.top_iterations.keys())
        while remaining_branches != 0:
            iteration = epoch.top_iterations[remaining_iterations.pop()]
            old_aa, new_aa = self.__get_random_aa__(iteration)
            # Build the mutation object for this iteration.
            # mut_chainID = iteration.chainIDs[self.idx_chain]
            # mut_resSeq = iteration.resSeqs[self.idx_chain][self.idx_residue]
            mutation = Mutation(
                chainID=mut_chainID,
                resSeq=mut_resSeq,
                old_aa=old_aa,
                new_aa=new_aa,
                chainID_idx=self.idx_chain,
                resSeq_idx=self.idx_residue,
            )

            self.mutations[iteration.iter_name].append(mutation)
            if len(remaining_iterations) == 0:
                # If all iterations have been mutated at least once and we still
                # have branches to generate, restart `remaining_iterations`.
                remaining_iterations = set(epoch.top_iterations.keys())
            remaining_branches -= 1

    def __get_random_pos__(self, epoch: Epoch) -> Tuple[str, int]:
        # Get an iteration to read the chainIDs and the resSeqs.
        any_iteration = next(iter(epoch.top_iterations.values()))

        # Choose the position to mutate. This will be the same for all iterations.
        max_tries: int = sum([len(resSeq) for resSeq in any_iteration.resSeqs]) * 5
        for i in range(max_tries):
            n_chains = len(any_iteration.chainIDs)
            idx_chain = choice(range(0, n_chains))
            n_residues = len(any_iteration.resSeqs[self.idx_chain])
            idx_residue = choice(range(0, n_residues))

            mut_resSeq = any_iteration.resSeqs[idx_chain][idx_residue]
            if mut_resSeq not in self.excluded_pos:
                self.idx_chain = idx_chain
                self.idx_residue = idx_residue
                mut_chainID = any_iteration.chainIDs[self.idx_chain]
                return mut_chainID, mut_resSeq

        raise RuntimeError(
            f"Can't generate novel position after {max_tries} tries. This shouldn't happen."
            f"Excluded positions: {self.excluded_pos}."
        )

    def __pop_random_aa__(self) -> str:
        cat_idx = choice(tuple(self.remaining_categories))
        self.remaining_categories.difference_update({cat_idx})

        shuffled_aas = sample(self.CAT_AAS[cat_idx], len(self.CAT_AAS[cat_idx]))
        for new_aa in shuffled_aas:
            if new_aa not in self.excluded_aas:
                self.excluded_aas.add(new_aa)
                return new_aa
        raise RuntimeError("Can't generate novel AA. This shouldn't happen.")

    def __get_random_aa__(self, iter: Iteration) -> Tuple[str, str]:
        old_aa = iter.resnames[self.idx_chain][self.idx_residue]
        self.excluded_aas.add(old_aa)
        new_aa = self.__pop_random_aa__()

        if len(self.remaining_categories) == 0:
            # All categories have already been chosen from `N` times.
            # Allow all of them again for the `N+1` iteration, except those
            # that are exhausted already.
            self.remaining_categories = set(range(self.N_CAT))
            for i in range(self.N_CAT):
                if set(self.CAT_AAS[i]).issubset(self.excluded_aas):
                    # All AAs from thise category have already been chosen
                    self.remaining_categories.difference_update({i})
        return old_aa, new_aa

    # def __generate_mutation__(iteration: Iteration) -> Mutation:
    #     pass

    def __getitem__(self, key: str) -> Sequence[Mutation]:
        return self.mutations[key]

    def __iter__(self) -> Iterator:
        return self.mutations.__iter__()

    def __contains__(self, value: Iteration) -> bool:
        return self.mutations.__contains__(value)

    def __len__(self) -> int:
        return len(self.mutations)

    def __items__(self) -> ItemsView:
        return self.mutations.items()
