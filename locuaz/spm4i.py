from collections.abc import Iterable, ItemsView
from typing import List, Tuple, Set, Iterator
from random import choice
from collections import defaultdict
from logging import Logger
from deprecated.sphinx import deprecated

from locuaz.projectutils import Branch, Epoch
from locuaz.mutation import Mutation
from locuaz.interface import get_interfacing_residues
from locuaz.spm4 import SPM4

@deprecated(version="0.7.0",
            reason="Mutation Generators are replaced by Mutation Creators.")
class SPM4i(SPM4):
    """
    overrides ``__generate_position__()``. It uses freesasa to get the interface to prevent mutations on positions
    that are not in contact.
    """
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
        self.probe_radius = probe_radius

        mut_idx_chain, mut_chainID, mut_idx_residue, mut_resSeq = self.__generate_position__(epoch, use_tleap, logger)
        self.__fill_mutations__(epoch, branches, mut_idx_chain, mut_chainID, mut_idx_residue, mut_resSeq)

    def __generate_position__(self, epoch: Epoch, use_tleap: bool = False,
                              logger: Logger = None) -> Tuple[int, str, int, int]:
        # Get a branch to read the chainIDs and the resSeqs.
        try:
            any_branch = next(iter(epoch.top_branches.values()))
        except Exception:
            raise RuntimeError(
                f"No available branches on Epoch {epoch.id}. "
                "It's likely that all of them failed during MD."
            )

        # Get the resSeq of all the residues in the interface and extend it
        interface_resSeq = get_interfacing_residues(any_branch.complex.pdb, any_branch.chainIDs, self.probe_radius,
                                                    use_tleap)

        # Now, filter the mutating resSeqs to keep only the residues that lie on the interface.
        candidates_resSeq: List[int] = []
        for cdr in any_branch.resSeqs:
            for resSeq in cdr:
                if resSeq in interface_resSeq and resSeq not in self.excluded_pos:
                    candidates_resSeq.append(resSeq)
        if len(candidates_resSeq) == 0:
            raise RuntimeError(f"No CDR residue on the interface that isn't excluded.")

        logger.info(f"Generating mutations with: {self}.\n"
                    f"resSeq at the interface: {interface_resSeq}.\n"
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
        return super().__pop_random_aa__()

    def __get_random_aa__(self, branch: Branch, idx_chain: int, idx_residue: int) -> Tuple[str, str]:
        return super().__get_random_aa__(branch, idx_chain, idx_residue)

    def __fill_mutations__(self, epoch: Epoch, branches: int, mut_idx_chain: int, mut_chainID: str,
                           mut_idx_residue: int, mut_resSeq: int):
        return super().__fill_mutations__(epoch, branches, mut_idx_chain, mut_chainID, mut_idx_residue, mut_resSeq)

    def __getitem__(self, key: str) -> Iterable[Mutation]:
        return super().__getitem__(key)

    def __iter__(self) -> Iterator:
        return super().__iter__()

    def __contains__(self, value: Branch) -> bool:
        return super().__contains__(value)

    def __len__(self) -> int:
        return super().__len__()

    def __items__(self) -> ItemsView:
        return super().__items__()

    def __str__(self) -> str:
        return "Single Point Mutation from 4 amino acid groups and interface detection "\
               f"with a probe radius of {self.probe_radius}"
