from collections.abc import Iterable, ItemsView
from typing import List, Tuple, Set, Iterator
from random import choice
from collections import defaultdict
from logging import Logger

from projectutils import Iteration, Epoch
from mutator import Mutation
from interface import get_interfacing_residues
from spm4i import SPM4i


class SPM4ie(SPM4i):

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
        self.__fill_mutations__(epoch, branches, mut_idx_chain, mut_chainID, mut_idx_residue, mut_resSeq)

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

        # Get the resSeq of all the residues in the interface and extend it
        interface_resSeq = get_interfacing_residues(any_iteration.complex.pdb, any_iteration.chainIDs, use_tleap)
        extended_interface_resSeq: Set[int] = set()
        for resSeq in interface_resSeq:
            extended_interface_resSeq.add(resSeq-1)
            extended_interface_resSeq.add(resSeq)
            extended_interface_resSeq.add(resSeq+1)


        # Now, filter the mutating resSeqs to keep only the residues that lie on the interface.
        mutating_interface_resSeq: List[int] = []
        for cdr in any_iteration.resSeqs:
            for resSeq in cdr:
                if resSeq in extended_interface_resSeq and resSeq not in self.excluded_pos:
                    mutating_interface_resSeq.append(resSeq)
        if len(mutating_interface_resSeq) == 0:
            raise RuntimeError(f"No CDR residue on the interface that isn't excluded.")

        logger.info(f"Generating mutations with: {self}.\n"
                    f"resSeq at the extended interface: {extended_interface_resSeq}.\n"
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
        return super().__pop_random_aa__()

    def __get_random_aa__(self, iteration: Iteration, idx_chain: int, idx_residue: int) -> Tuple[str, str]:
        return super().__get_random_aa__(iteration, idx_chain, idx_residue)

    def __fill_mutations__(self, epoch: Epoch, branches: int, mut_idx_chain: int, mut_chainID: str,
                           mut_idx_residue: int, mut_resSeq: int):
        return super().__fill_mutations__(epoch, branches, mut_idx_chain, mut_chainID, mut_idx_residue, mut_resSeq)

    def __getitem__(self, key: str) -> Iterable[Mutation]:
        return super().__getitem__(key)

    def __iter__(self) -> Iterator:
        return super().__iter__()

    def __contains__(self, value: Iteration) -> bool:
        return super().__contains__(value)

    def __len__(self) -> int:
        return super().__len__()

    def __items__(self) -> ItemsView:
        return super().__items__()

    def __str__(self) -> str:
        return "Single Point Mutation from 4 amino acid groups and use of the interface"