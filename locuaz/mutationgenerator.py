from abc import ABC, abstractmethod
from email import iterators
from typing import List, Tuple, Dict, Set, Iterator
from random import choice, sample
from collections import defaultdict
from collections.abc import ItemsView, Mapping

from projectutils import Iteration, Epoch
from mutator import Mutation

# fmt: off
# AA3_LIST = ("Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "His", "Ile",
#     "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp",  "Tyr", "Val")
# AA1_LIST = ('A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
#     'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
# # Cystein is not included.
# NEG_AA1_LIST = ('D', 'E', 'S', 'T')
# POS_AA1_LIST = ('R', 'N',  'Q', 'H', 'K')
# PHO_AA1_LIST = ('A', 'G', 'I', 'M', 'L', 'V')
# RIN_AA1_LIST = ('P', 'F', 'W', 'Y')
# CAT_AA1_LIST = (NEG_AA1_LIST, POS_AA1_LIST, PHO_AA1_LIST, RIN_AA1_LIST)
# fmt: on


class AbstractMutationGenerator(ABC):
    @abstractmethod
    def __init__(self, epoch: Epoch) -> None:
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
    def __init__(self, epoch: Epoch) -> None:
        pass

    def __getitem__(self, key: Iteration) -> Mutation:
        pass

    def __iter__(self) -> Iterator:
        pass

    def __contains__(self, value: Iteration) -> bool:
        pass


class SPM_4(AbstractMutationGenerator, Mapping):
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
    idx_chain: int
    idx_residue: int
    # remaining_branches: int
    # remaining_categories: set(CAT_AAS)

    # new_aas: Set
    # iterations: List[Iteration]

    mutations: Dict[str, List[Mutation]]

    def __get_random_aa__(self) -> str:
        cat_idx = choice(tuple(self.remaining_categories))
        self.remaining_categories.difference_update({cat_idx})

        shuffled_aas = sample(self.CAT_AAS[cat_idx], len(self.CAT_AAS[cat_idx]))
        for new_aa in shuffled_aas:
            if new_aa not in self.excluded_aas:
                self.excluded_aas.add(new_aa)
                return new_aa
        raise RuntimeError("Can't generate new binder. This shouldn't happen.")

    def __pop_random_aa__(self, iter: Iteration) -> Tuple[str, str]:
        old_aa = iter.resnames[self.idx_chain][self.idx_residue]
        self.excluded_aas.add(old_aa)
        new_aa = self.__get_random_aa__()

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

    def __generate_mutation__(iteration: Iteration) -> Mutation:
        pass

    def __init__(self, epoch: Epoch, max_branches: int) -> None:
        self.remaining_categories = set(range(0, self.N_CAT))
        self.mutations = defaultdict(list)
        self.excluded_aas = set()
        iteration = epoch[epoch.top_iterations[0]]

        # First, choose the chain and position to mutate. These will be
        # the same for all iterations.
        n_chains = len(iteration.chainIDs)
        self.idx_chain = choice(range(0, n_chains))
        n_residues = len(iteration.resSeqs[self.idx_chain])
        self.idx_residue = choice(range(0, n_residues))

        # Now, generate up to `max_branches` mutations
        remaining_branches = max_branches
        remaining_iterations = set(epoch.top_iterations)
        while remaining_branches != 0:
            iteration = epoch[remaining_iterations.pop()]
            old_aa, new_aa = self.__pop_random_aa__(iteration)
            # Build the mutation object for this iteration.
            mut_chainID = iteration.chainIDs[self.idx_chain]
            mut_resSeq = iteration.resSeqs[self.idx_chain][self.idx_residue]
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
                remaining_iterations = set(epoch.top_iterations)
            remaining_branches -= 1

    def __getitem__(self, key: Iteration) -> Mutation:
        return self.mutations[key]

    def __iter__(self) -> Iterator:
        return self.mutations.__iter__()

    def __contains__(self, value: Iteration) -> bool:
        return self.mutations.__contains__(value)

    def __len__(self) -> int:
        return len(self.mutations)

    def __items__(self) -> ItemsView:
        return self.mutations.items()


# def generate_mutations(iteration: Iteration, *, max_branches: int) -> List[Mutation]:

#     # First, choose the position to mutate:
#     n_chains = len(iteration.chainIDs)
#     idx_chain = choice(range(0, n_chains))
#     n_residues = len(iteration.resSeqs[idx_chain])
#     idx_residue = choice(range(0, n_residues))
#     # Then, save the current AA
#     old_aa = iteration.resnames[idx_chain][idx_residue]

#     n_cat = len(CAT_AA1_LIST)
#     categories = set(range(n_cat))
#     new_aas = {old_aa}
#     while max_branches != 0:
#         cat_idx = choice(tuple(categories))
#         categories.difference_update({cat_idx})

#         shuffled_aas = sample(CAT_AA1_LIST[cat_idx], len(CAT_AA1_LIST[cat_idx]))
#         for new_aa in shuffled_aas:
#             if new_aa not in new_aas:
#                 new_aas.add(new_aa)
#                 max_branches -= 1
#                 break
#         else:
#             raise RuntimeError(
#                 "Can't generate new binder. This is a logic error. "
#                 "This shouldn't happen."
#             )

#         if len(categories) == 0:
#             # All categories have already been chosen from `N` times.
#             # Allow all of them again for the `N+1` iteration, except those
#             # that are exhausted already
#             categories = set(range(n_cat))
#             for i in range(n_cat):
#                 if set(CAT_AA1_LIST[i]).issubset(new_aas):
#                     # All AAs from thise category have already been chosen
#                     categories.difference_update({i})

#     # Finally, build the list of mutation objects
#     new_aas.difference_update({old_aa})
#     mut_chainID = iteration.chainIDs[idx_chain]
#     mut_resSeq = iteration.resSeqs[idx_chain][idx_residue]
#     mutations = [
#         Mutation(
#             chainID=mut_chainID,
#             resSeq=mut_resSeq,
#             old_aa=old_aa,
#             new_aa=aa,
#             chainID_idx=idx_chain,
#             resSeq_idx=idx_residue,
#         )
#         for aa in new_aas
#     ]

#     return mutations
