from collections.abc import Iterable, ItemsView
from typing import List, Tuple, Set, Iterator, Dict, Optional
from collections import defaultdict
from logging import Logger
from pathlib import Path
import csv
from deprecated.sphinx import deprecated

from locuaz.projectutils import Branch, Epoch
from locuaz.mutation import Mutation
from locuaz.interface import get_interfacing_residues
from locuaz.spm4 import SPM4

@deprecated(version="0.7.0",
            reason="Mutation Generators are replaced by Mutation Creators.")
class SPM4gmxmmpbsa(SPM4):
    """
        overrides ``__generate_position__()``. It uses freesasa to get the interface to prevent mutations on positions
        that are not in contact and also sorts residues according to the binding ΔG, reported by the gmxmmpbsa
        scoring function
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
                              logger: Optional[Logger] = None) -> Tuple[int, str, int, int]:
        # Get a branch to read the chainIDs and the resSeqs.
        try:
            any_branch = next(iter(epoch.top_branches.values()))
        except Exception:
            raise RuntimeError(f"No available branches on Epoch {epoch.id}. "
                               "It's likely that all of them failed during MD.")

        # Get the resSeq of all the residues in the interface and extend it
        interface_resSeq = get_interfacing_residues(any_branch.complex.pdb, any_branch.chainIDs,
                                                    self.probe_radius,
                                                    use_tleap)

        # Now, filter the mutating resSeqs to keep only the residues that lie on the interface.
        precandidates_resSeq: Set[int] = set()
        for cdr in any_branch.resSeqs:
            for resSeq in cdr:
                if resSeq in interface_resSeq and resSeq not in self.excluded_pos:
                    precandidates_resSeq.add(resSeq)
        if len(precandidates_resSeq) == 0:
            raise RuntimeError(f"No CDR residue on the interface that isn't excluded.")

        candidates_resSeq = self.__sort_by_deltas__(precandidates_resSeq, epoch, logger)

        logger.info(f"Generating mutations with: {self}.\n"
                    f"resSeq at the interface: {interface_resSeq}.\n"
                    f"'mutating_resSeq': {any_branch.resSeqs}.\n"
                    f"excluded resSeq: {self.excluded_pos}.\n"
                    f"'mutating_resSeq' that may be mutated, inversely ordered by ΔG: {candidates_resSeq}.")

        # Choose the position to mutate. This will be the same for all branches.
        mut_resSeq = self.__choose_resSeq__(candidates_resSeq)
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

    def __sort_by_deltas__(self, resSeqs: Set[int], epoch: Epoch, logger: Optional[Logger] = None) -> List[int]:
        """
        __sort_by_deltas__() accumulates the ΔG of each CDR residue given by gmxmmpbsa idecomp and sorts the input
        candidate resSeqs so the least interacting are at the top.
        Args:
            resSeqs (Set[int]): set of candidate resSeqs to be mutated
            epoch (Epoch): old epoch with top_branches

        Returns:
            List[int]: candidate resSeqs sorted by increasing ΔG
        """
        deltas = defaultdict(int)
        for branch_name, branch in epoch.top_branches.items():
            decomp_mmpbsa = Path(branch.score_dir, "gmxmmpbsa", "decomp_gmxmmpbsa.csv")
            assert decomp_mmpbsa.is_file(), f"No 'decomp_gmxmmpbsa.csv' in {Path(branch.score_dir, 'gmxmmpbsa')}. " \
                                            f"Cannot generate mutation. "
            delta_iter = self.__get_deltas__(decomp_mmpbsa)
            for resSeq, delta_G in delta_iter.items():
                deltas[resSeq] += delta_G
        sorted_deltas: List[Tuple[int, float]] = sorted(deltas.items(), key=lambda item: item[1], reverse=True)
        self.__log_deltas__(sorted_deltas, logger)

        sorted_resSeq: List[int] = []
        for resSeq, delta_g in sorted_deltas:
            if resSeq in resSeqs:
                resSeqs.remove(resSeq)
                sorted_resSeq.append(resSeq)
        # `resSeqs` should be empty by now, but just in case, add the remaining at the end.
        if len(resSeqs) > 0:
            for resSeq in resSeqs:
                sorted_resSeq.append(resSeq)
        return sorted_resSeq

    @staticmethod
    def __get_deltas__(csv_path: Path) -> Dict[int, float]:
        """
        __get_deltas__(): reads a .csv file from gmxmmpbsa with the ΔG residue decompostion from gmxmmpbsa
        Args:
            csv_path (Path): gmxmmpbsa csv decomposition file

        Returns:
            defaultdict[int, float]: resSeq and its ΔG
        """
        with open(Path(csv_path), 'r') as csv_file:
            text = csv.reader(csv_file)
            for line in text:
                if line == ["DELTAS:"]:
                    next(text)
                    next(text)
                    break
                else:
                    continue
            deltas = defaultdict(int)
            for linea in text:
                if linea[0] == '\n':
                    break
                rl, chainID, resname, resSeq = linea[1].split(':')
                if chainID == 'B':
                    # deltas[(resSeq, chainID, resname)] = float(linea[-1])
                    # noinspection PyTypeChecker
                    deltas[int(resSeq)] = float(linea[-1])
        return deltas

    @staticmethod
    def __log_deltas__(deltas: List[Tuple[int, float]], logger: Optional[Logger] = None) -> None:
        if logger:
            for i in range(len(deltas)):
                deltas[i] = (deltas[i][0], round(deltas[i][1], 3))
            logger.info(f"(resSeq, ΔG) tuples: {deltas}.")

    @staticmethod
    def __choose_resSeq__(candidates_resSeq: List[int]) -> int:
        return candidates_resSeq[0]

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
        return "Single Point Mutation from 4 amino acid groups, ordered by gmxmmpbsa, and interface detection " \
               f"with a probe radius of {self.probe_radius}"
