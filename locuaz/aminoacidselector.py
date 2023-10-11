from collections import defaultdict
from logging import Logger
from warnings import warn
from random import choice, choices
from typing import List, Dict, Set, Any, Final, Tuple, Union, Optional

from locuaz.mutation import Mutation
from locuaz.projectutils import Branch
from locuaz.siteselector import Site


class AminoAcidSelector:
    use_bins: bool
    bins: Dict[int, Set[str]]
    aa_to_bin: Dict[str, Set[int]] = {"": set()}
    N_BINS: int
    N_SITES: int
    bins_criteria: str
    uniform: Dict[str, float] = {
        'A': 0.05263, 'R': 0.05263, 'N': 0.05263, 'D': 0.05263, 'E': 0.05263,
        'Q': 0.05263, 'G': 0.05263, 'H': 0.05263, 'I': 0.05263, 'L': 0.05263,
        'K': 0.05263, 'M': 0.05263, 'F': 0.05263, 'P': 0.05263, 'S': 0.05263,
        'T': 0.05263, 'W': 0.05263, 'Y': 0.05263, 'V': 0.05263, 'C': 0.00000}
    # Currently not used
    reis_barletta_full: Final[Dict[str, float]] = {
        'A': 0.03353, 'R': 0.05076, 'N': 0.06358, 'D': 0.06073, 'E': 0.02587,
        'Q': 0.01545, 'G': 0.07939, 'H': 0.0219, 'I': 0.03489, 'L': 0.0385,
        'K': 0.0185, 'M': 0.00678, 'F': 0.04677, 'P': 0.0234, 'S': 0.12669,
        'T': 0.06186, 'W': 0.05624, 'Y': 0.20396, 'V': 0.0312, 'C': 0.00000}
    reis_barletta_cdr: Final[Dict[str, float]] = {
        'A': 0.03579, 'R': 0.04462, 'N': 0.06788, 'D': 0.06806, 'E': 0.02621,
        'Q': 0.01268, 'G': 0.09396, 'H': 0.02256, 'I': 0.02952, 'L': 0.0387,
        'K': 0.01514, 'M': 0.0076, 'F': 0.05339, 'P': 0.02284, 'S': 0.12746,
        'T': 0.06049, 'W': 0.0493, 'Y': 0.19701, 'V': 0.02678, 'C': 0.00000}
    aa_distribution: Dict[str, float] = uniform

    excluded_aas: Set[str]
    full_bin: Final[Tuple[str]] = ('A', 'R', 'N', 'D', 'E', 'Q', 'G', 'H', 'I',
                                   'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y',
                                   'V')

    def __init__(self, excluded_aas: Set[str], creation_config: Dict[str, Any]) -> None:
        self.bins = {}
        for i, each_bin in enumerate(creation_config["aa_bins"]):
            # Store bins
            self.bins[i] = set(each_bin)
            # Store in which bin is each amino acid
            for aa in each_bin:
                self.aa_to_bin[aa] = set([i])
        self.N_BINS = len(self.bins)

        self.bins_criteria = creation_config["aa_bins_criteria"]
        # TODO: move this up to cli.py and make distributions csv datafiles.
        if creation_config["aa_probability"] == "ReisBarletta":
            self.aa_distribution = self.reis_barletta_cdr
        elif creation_config["aa_probability"] == "custom":
            self.aa_distribution = creation_config["custom"]
        else:
            # It's uniform
            pass

        self.excluded_aas = excluded_aas

    def __call__(self,
                 top_branches: Dict[str, Branch],
                 branches: int,
                 sites: List[Site],
                 log: Optional[Logger]) -> Dict[str, List[Mutation]]:
        """
        Generate mutations according to user input (bins, probability of
        choosing each amino acid, etc...).
        Parameters
        ----------
        top_branches : Dict[str, Branch]
            top branches from epoch `i`.
        branches : int
            number of new branches to generate in total.
        sites : List[Site]
            sites where each branch will be mutated.
        Returns
        -------
        Dict[str, List[Mutation]]
            dictionary with branch names from epoch `i` and the mutations that
            will be applied on them to generate the branches of epoch `i+1`.
        """
        mutations = defaultdict(list)
        new_branches = branches
        remaining_branches = set(top_branches.keys())
        remaining_bins = set()

        while new_branches != 0:
            # prevent nbr of sites being higher than nbr of new branches
            # And probably larger that nbr of new branches * nbr of bins.
            for site in sites:
                branch = top_branches[remaining_branches.pop()]
                old_aa = branch.resnames[site.idx_chain][site.idx_residue]
                try:
                    new_aa = self.__select_aa__(old_aa, remaining_bins)
                except AssertionError as e:
                    log.warning(f"{e}.\nGenerated {branches-new_branches} branches.")
                    break
                # Add this mutation to this branch.
                mutations[branch.branch_name].append(Mutation.from_site(
                    site, old_aa=old_aa, new_aa=new_aa)
                )
                if len(remaining_branches) == 0:
                    # If all branches have been mutated at least once, and we still
                    # have branches to generate, restart `remaining_branches`.
                    remaining_branches = set(top_branches.keys())
                new_branches -= 1
        return mutations

    def __select_aa__(self, old_aa: str, remaining_bins: Set[int]) -> str:
        self.excluded_aas.add(old_aa)

        if len(remaining_bins) == 0:
            remaining_bins = self.__reset_bins__(old_aa)
        # Choose bin
        bin_idx = choice(tuple(remaining_bins))
        remaining_bins.difference_update([bin_idx])

        # Get a new amino acid.
        try:
            aa = self.__pop_random_aa__(bin_idx)
            return aa
        except AssertionError as e:
            raise e

    def __reset_bins__(self, old_aa: str) -> Set[int]:
        """
        If the bins_criteria is set to "exclusive":
            when all bins (amino acid categories) have already been chosen from,
            allow all of them again for the next branches.
            If all the amino acids form a certain bin have been chosen already,
            then that bin will not be reset. Also discard the bin corresponding
            to the ``old_aa``.
        If the bins_criteria is set to "within":
            Just return the bin corresponding to ``old_aa``.
        Parameters
        ----------
        old_aa : str
            Current AA in the chosen position.
        Returns
        -------
        remaining_bins: Set[int]
            reset bins
        """
        if self.bins_criteria == "exclusive":
            remaining_bins = set(range(0, self.N_BINS))
            remaining_bins.difference_update(self.aa_to_bin[old_aa])
            for i in range(self.N_BINS):
                if self.bins[i].issubset(self.excluded_aas):
                    # All AAs from this category have already been chosen
                    remaining_bins.difference_update({i})
        elif self.bins_criteria == "within":
            # Make sure to deep copy the set that's inside the dict
            remaining_bins = set(self.aa_to_bin[old_aa])
        else:
            raise ValueError(f"This shouldn't happen. {self.bins_criteria=}")
        return remaining_bins

    def __pop_random_aa__(self, bin_idx: int) -> str:
        current_bin = self.bins[bin_idx]
        while True:
            norm_bin = self.__get_norm_bin__(self.aa_distribution, current_bin)
            new_aa = choices(list(norm_bin.keys()), list(norm_bin.values()))[0]
            if new_aa not in self.excluded_aas:
                self.excluded_aas.add(new_aa)
                return new_aa
            else:
                current_bin.difference_update(set(new_aa))
                assert len(current_bin) > 0, "Could not choose aa from " \
                                             f"{self.bins[bin_idx]}. All AAs  " \
                                             f"are excluded. {self.excluded_aas=}"

    def __get_norm_bin__(self,
                         distro: Dict[Union[str, int], float],
                         aa_bin: Set[str]) -> Dict[Union[str, int], float]:
        return self.__norm_distribution__(
            {aa: prob for aa, prob in distro.items() if aa in aa_bin})

    @staticmethod
    def __norm_distribution__(distro: Dict[Union[str, int], float]
                              ) -> Dict[Union[str, int], float]:
        tot = sum(distro.values())
        return {aa: prob / tot for aa, prob in distro.items()}
