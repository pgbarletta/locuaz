import csv
from collections import defaultdict
from itertools import chain
from logging import Logger
from pathlib import Path
from random import shuffle
from typing import List, Dict, Set, Final, Any, Optional, Tuple

from locuaz.interface import get_interfacing_residues
from locuaz.projectutils import Branch
from locuaz.mutation import Site

class SiteSelector:
    """
    SiteSelector chooses the positions to mutate by:

    * taking the mutating resSeq.
    * excluding those positions that are in ``excluded_resSeqs``
    * optionally excluding those positions that are not in the target-binder
      interface
    * randomly sorting the remaining position or sorting them according to
      their contribution to the free energy of binding according to the
      output from ``gmxmmpbsa`` so the less contributing are picked first.

    Then it will return a list of Sites as long as requested.

    Parameters
    ----------
    all_resSeqs : List[List[int]]
        positions that may be mutated. They may be grouped into different lists,
        without this having any effect.
    chainIDs : List[str]
        chainIDs of the target chains that may be mutated
    excluded_resSeqs : Set[int]
        resSeqs of the excluded positions, probably because they've recently
        been mutated
    creation_config : Dict[str, Any]
        creation section from the input config file
    amber_numbering : bool
        when using Tleap, the resSeq numbering scheme is continuous as opposed
        to the strided scheme from GROMACs where each chain begins with resSeq 1
    """
    N_SITES: int
    all_resSeqs: List[List[int]]
    chainIDs: List[str]
    excluded_resSeqs: Set[int]
    probability: str
    use_interface: bool
    probe_radius: Optional[float]
    amber_numbering: bool

    idx_chain: int
    idx_residue: int

    def __init__(self, all_resSeqs: List[List[int]], chainIDs: List[str],
                 excluded_resSeqs: Set[int],
                 creation_config: Dict[str, Any], *, amber_numbering: bool = False) -> None:
        self.all_resSeqs = all_resSeqs
        self.chainIDs = chainIDs
        self.excluded_resSeqs = excluded_resSeqs

        self.N_SITES = creation_config["sites"]
        self.probability = creation_config["sites_probability"]
        self.use_interface = creation_config["sites_interfacing"]
        self.probe_radius = creation_config["sites_interfacing_probe_radius"]

        self.amber_numbering = amber_numbering

    def __call__(self, top_branches: Dict[str, Branch], logger: Optional[Logger] = None) -> List[Site]:

        candidates_resSeq = self.__filter_candidates__(top_branches, logger)
        sites: List[Site] = []
        for i in range(self.N_SITES):
            # Choose the position to mutate.
            mut_resSeq = candidates_resSeq[i]
            # Now, get the remaining details associated to `mut_resSeq`,
            # including the chain from where it came from
            for j, resSeqs in enumerate(self.all_resSeqs):
                if mut_resSeq in resSeqs:
                    sites.append(Site(idx_chain=j,
                                      chainID=self.chainIDs[j],
                                      idx_residue=resSeqs.index(mut_resSeq),
                                      resSeq=mut_resSeq))
                    break
            else:
                raise RuntimeError(f"Tried to mutate resSeq {mut_resSeq} but "
                                   "could not find it in the mutating resSeqs: "
                                   f"{self.all_resSeqs}\n"
                                   "This should not happen.")
        return sites

    def __filter_candidates__(self, top_branches: Dict[str, Branch],
                              logger: Optional[Logger] = None) -> List[int]:
        if self.use_interface:
            # Get the resSeq of all the residues in the interfaces of all branches.
            interface_resSeq: Set[int] = set()
            for branch in top_branches.values():
                interface_resSeq |= get_interfacing_residues(branch.complex.pdb,
                                                             self.chainIDs,
                                                             self.probe_radius,
                                                             self.amber_numbering)
                assert len(interface_resSeq), \
                    f"{branch} has no interfacing residues. It cannot be  top branch."
            allowed_resSeq = interface_resSeq
        else:
            # Pretend as if everything is on the interface.
            allowed_resSeq = {resSeq for cdr in self.all_resSeqs for resSeq in cdr}

        # Remove the excluded ones.
        allowed_resSeq.difference_update(self.excluded_resSeqs)
        # Now, filter the mutating resSeqs.
        candidates_resSeq = set(chain.from_iterable(self.all_resSeqs)).intersection(allowed_resSeq)

        if len(candidates_resSeq) == 0:
            err_msg = f"No CDR residue on the interface that isn't excluded.\n" \
                      f"Excluded positions: {self.excluded_resSeqs}"
            if self.use_interface:
                err_msg += f"\nInterfacing positions: {interface_resSeq}"
            raise RuntimeError(err_msg)

        if self.probability == "mmpbsa":
            sorted_candidates_resSeq: List[int] = self.__sort_by_deltas__(
                candidates_resSeq, top_branches, logger)
        else:
            sorted_candidates_resSeq: List[int] = list(candidates_resSeq)
            shuffle(sorted_candidates_resSeq)

        return sorted_candidates_resSeq

    def __sort_by_deltas__(self, resSeqs: Set[int], top_branches: Dict[str, Branch],
                           logger: Optional[Logger] = None) -> List[int]:
        """
        accumulates the ΔG of each CDR residue given by gmxmmpbsa idecomp and sorts the input
        candidate resSeqs so the least interacting are at the top.
        Args:
            resSeqs (Set[int]): set of candidate resSeqs to be mutated
            top_branches (Dict[str, Branch]): top_branches from the previous epoch

        Returns:
            List[int]: candidate resSeqs sorted by increasing ΔG
        """
        deltas = defaultdict(int)
        for branch_name, branch in top_branches.items():
            decomp_mmpbsa = Path(branch.score_dir, "gmxmmpbsa", "decomp_gmxmmpbsa.csv")
            assert decomp_mmpbsa.is_file(), "No 'decomp_gmxmmpbsa.csv' in " \
                                            f"{Path(branch.score_dir, 'gmxmmpbsa')}. " \
                                            "Cannot select mutation site. "
            residue_deltas = self.__get_deltas__(decomp_mmpbsa)
            for resSeq, delta_G in residue_deltas.items():
                deltas[resSeq] += delta_G
        sorted_deltas: List[Tuple[int, float]] = sorted(
            deltas.items(), key=lambda item: item[1], reverse=True)
        self.__log_deltas__(sorted_deltas, logger)

        sorted_resSeq: List[int] = []
        for resSeq, delta_g in sorted_deltas:
            if resSeq in resSeqs:
                resSeqs.remove(resSeq)
                sorted_resSeq.append(resSeq)
        # `resSeqs` should be empty by now, but just in case, add the remaining at the end.
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
