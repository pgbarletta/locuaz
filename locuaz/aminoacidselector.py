from collections import defaultdict
from copy import deepcopy
from logging import Logger
from random import choices
from typing import List, Dict, Set, Any, Final, Tuple, Union

from locuaz.mutation import Mutation
from locuaz.projectutils import Branch
from locuaz.siteselector import Site
from pools import BinPool, InfinitePool
from primitives import GenerationError


class AminoAcidMemory:
    # Branch names + site are the keys in these 2 dicts.
    # We use them to better decide the next mutation on each branch.
    branch_bins: Dict[Tuple[str, Site], Dict[int, Set[str]]]
    branch_bins_indices: Dict[Tuple[str, Site], BinPool]
    epoch_bins_indices: BinPool

    def __init__(self, top_branches: Dict[str, Branch],
                 bins: Dict[int, Set[str]],
                 sites: List[Site]) -> None:

        self.nbins = len(bins)
        self.epoch_bins_indices = BinPool.from_size(self.nbins)

        self.branch_bins = {}
        self.branch_bins_indices = {}
        for site in sites:
            for branch_name in top_branches.keys():
                self.branch_bins[(branch_name, site)] = deepcopy(bins)
                self.branch_bins_indices[(branch_name, site)] = BinPool.from_size(self.nbins)

    def remove_from_memory(self, branch_name: str, site: Site, bin_idx: int, aa: str) -> None:
        self.branch_bins_indices[(branch_name, site)].discard(bin_idx)
        self.epoch_bins_indices.discard(bin_idx)
        self.branch_bins[(branch_name, site)][bin_idx].discard(aa)

    def get_bin(self, branch_name: str, site: Site,
                old_aa: str, excluded_bins: Set[int]) -> int:

        while True:
            bin_idx = self.__attempt_novel_bin_idx__(branch_name, site, excluded_bins)
            bin_aas = self.branch_bins[(branch_name, site)][bin_idx]
            bin_aas.discard(old_aa)
            if len(bin_aas) > 0:
                return bin_idx
            excluded_bins.add(bin_idx)
            if len(excluded_bins) == self.nbins:
                raise GenerationError(f"Cannot mutate branch {branch_name}. "
                                      "All bins have been excluded.")

    def __attempt_novel_bin_idx__(self,
                                  branch_name: str,
                                  site: Site,
                                  excluded_bins: Set[int]) -> int:
        """
        Tries to pick a random bin number that's available both in the branch
        and the epoch bin pools and isn't excluded.
        The pools are never empty, so if it can't find one, it's just that they
        have no common element that hasn't been excluded. In that case, it
        returns a bin index from the branch's bin pool.

        Parameters
        ----------
        branch_name : str
            branch that's being mutated
        site : Site
            mutation site
        excluded_bins : Set[int]


        Returns
        -------
        int
            index of a bin that's both available for the branch and for the
            epoch.
        Raises
        ------
        RunTimeError

        """
        common_bins: Set[int]
        branch_pool = self.branch_bins_indices[(branch_name, site)]
        branch_pool.difference_update(excluded_bins)
        self.epoch_bins_indices.difference_update(excluded_bins)

        common_bins = set(branch_pool & self.epoch_bins_indices)
        if len(common_bins) == 0:
            # `branch_name`'s bin pool and the epoch's bin pool have no common
            # elements that haven't been excluded.
            return branch_pool.pick()
        else:
            return BinPool(common_bins).pick()


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
    full_bin: Final[Tuple[str]] = ('A', 'R', 'N', 'D', 'E', 'Q', 'G', 'H', 'I',
                                   'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y',
                                   'V', 'C')

    aa_distribution: Dict[str, float] = uniform

    memory: AminoAcidMemory

    def __init__(self, creation_config: Dict[str, Any]) -> None:
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

    def __call__(self,
                 top_branches: Dict[str, Branch],
                 branches: int,
                 sites: List[Site],
                 logger: Logger) -> Dict[str, List[Mutation]]:
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
        old_branches_pool = InfinitePool(top_branches.keys())
        self.memory = AminoAcidMemory(top_branches, self.bins, sites)

        while new_branches != 0:
            branch = top_branches[old_branches_pool.pop()]
            for site in sites:
                old_aa = branch.resnames[site.idx_chain][site.idx_residue]
                try:
                    new_aa = self.__select_aa__(branch.branch_name, site, old_aa)
                    # Add this mutation to this branch.
                    mutations[branch.branch_name].append(Mutation.from_site(
                        site, old_aa=old_aa, new_aa=new_aa)
                    )
                except GenerationError as e:
                    logger.warning(f"{e}")

                new_branches -= 1
        return mutations

    def __select_aa__(self, branch_name: str, site: Site, old_aa: str) -> str:

        if self.bins_criteria == "exclude":
            excluded_bins: Set[int] = set(self.aa_to_bin[old_aa])
        else:
            # 'within'. All bins but the one that corresponds to the aa are excluded.
            excluded_bins: Set[int] = set(range(0, len(self.bins))) - self.aa_to_bin[old_aa]

        bin_idx = self.memory.get_bin(branch_name, site, old_aa, excluded_bins)
        # Normalize the probabilities of the AAs in the bin, so they add up to 1
        bin_aas = self.__get_norm_bin__(self.aa_distribution,
                                        self.memory.branch_bins[(branch_name, site)][bin_idx])
        new_aa = choices(list(bin_aas.keys()), list(bin_aas.values()))[0]
        self.memory.remove_from_memory(branch_name, site, bin_idx, new_aa)

        return new_aa

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
