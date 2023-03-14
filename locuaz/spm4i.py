from collections.abc import Iterable
from typing import List, Tuple, Dict, Set, Iterator
from random import choice, sample
from collections import defaultdict
from collections.abc import ItemsView
from logging import Logger

from projectutils import Iteration, Epoch
from mutator import Mutation
from interface import get_interfacing_residues
from abstractmutationgenerator import AbstractMutationGenerator


class SPM4i(AbstractMutationGenerator):
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
            branches: int,
            *,
            excluded_aas: Set[str],
            excluded_pos: Set[int],
            use_tleap: bool = False,
            logger: Logger = None
    ) -> None:
        self.remaining_categories = set(range(0, self.N_CAT))
        self.mutations = defaultdict(list)
        self.excluded_aas = excluded_aas
        self.excluded_pos = excluded_pos

        mut_idx_chain, mut_chainID, mut_idx_residue, mut_resSeq = self.__generate_position__(epoch, use_tleap, logger)
        # Now, generate up to `branches` mutations
        remaining_branches = branches
        remaining_iterations = set(epoch.top_iterations.keys())
        while remaining_branches != 0:
            iteration = epoch.top_iterations[remaining_iterations.pop()]
            old_aa, new_aa = self.__get_random_aa__(iteration, mut_idx_chain, mut_idx_residue)
            # Build the mutation object for this iteration.
            mutation = Mutation(
                chainID=mut_chainID,
                resSeq=mut_resSeq,
                old_aa=old_aa,
                new_aa=new_aa,
                chainID_idx=mut_idx_chain,
                resSeq_idx=mut_idx_residue,
            )

            self.mutations[iteration.iter_name].append(mutation)
            if len(remaining_iterations) == 0:
                # If all iterations have been mutated at least once, and we still
                # have branches to generate, restart `remaining_iterations`.
                remaining_iterations = set(epoch.top_iterations.keys())
            remaining_branches -= 1

    def __generate_position__(self, epoch: Epoch, use_tleap: bool = False,
                              logger: Logger = None) -> Tuple[int, str, int, int]:
        # Get an iteration to read the chainIDs and the resSeqs.
        try:
            any_iteration = next(iter(epoch.top_iterations.values()))
        except Exception:
            raise RuntimeError(
                f"No available iterations on Epoch {epoch.id}. "
                "It's likely that all of them failed during MD."
            )

        # Get the resSeq of all the residues in the interface
        interface_resSeq = get_interfacing_residues(any_iteration.complex.pdb, any_iteration.chainIDs, use_tleap)

        # Now, filter the mutating resSeqs to keep only the residues that lie on the interface.
        mutating_interface_resSeq: List[int] = []
        for cdr in any_iteration.resSeqs:
            for resSeq in cdr:
                if resSeq in interface_resSeq and resSeq not in self.excluded_pos:
                    mutating_interface_resSeq.append(resSeq)
        if len(mutating_interface_resSeq) == 0:
            raise RuntimeError(f"No CDR residue on the interface that isn't excluded.")

        logger.info(f"Generating mutations with: {self}.\n"
                    f"resSeq at the interface: {interface_resSeq}.\n"
                    f"'mutataing_resSeq': {any_iteration.resSeqs}.\n"
                    f"excluded resSeq: {self.excluded_pos}.\n"
                    f"'mutataing_resSeq' that may be mutated: {mutating_interface_resSeq}.")

        # Choose the position to mutate. This will be the same for all iterations.
        mut_resSeq = choice(mutating_interface_resSeq)
        # Now, get the remaining details associated to `mut_resSeq`, including the chain from where it came from
        for j, resSeqs in enumerate(any_iteration.resSeqs):
            if mut_resSeq in resSeqs:
                mut_idx_chain = j
                mut_chainID = any_iteration.chainIDs[mut_idx_chain]
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

    def __get_random_aa__(self, iteration: Iteration, idx_chain: int, idx_residue: int) -> Tuple[str, str]:
        old_aa = iteration.resnames[idx_chain][idx_residue]
        self.excluded_aas.add(old_aa)
        new_aa = self.__pop_random_aa__()

        if len(self.remaining_categories) == 0:
            # All categories have already been chosen from `N` times.
            # Allow all of them again for the `N+1` iteration, except those
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

    def __contains__(self, value: Iteration) -> bool:
        return self.mutations.__contains__(value)

    def __len__(self) -> int:
        return len(self.mutations)

    def __items__(self) -> ItemsView:
        return self.mutations.items()

    def __str__(self) -> str:
        return "Single Point Mutation from 4 amino acid groups and use of the interface."
