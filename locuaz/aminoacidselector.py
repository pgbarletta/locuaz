from collections import defaultdict
from copy import deepcopy
from itertools import chain
from logging import Logger
from pathlib import Path
from random import choices
from typing import List, Dict, Set, Any, Tuple, Union
import pandas as pd

from locuaz.mutation import Mutation
from locuaz.projectutils import Branch
from locuaz.siteselector import Site
from locuaz.pools import BinPool, InfinitePool
from locuaz.primitives import GenerationError


class AminoAcidMemory:
    """
    Class used to annotate the recently mutations created to avoid repeating
    them and to guide the next mutations for a better sampling of the solution
    space.

    Parameters
    ----------
    top_branches : Dict[str, Branch]
        top branches from the last epoch
    bins : Dict[int, Set[str]]
        Bins are numbered starting from 0 (the keys). Each bin (the values) is a
        Set of strings, each a one-letter coded amino acid.
    sites : List[Site]
        sites that will be mutated.
    Attributes
    ----------
    branch_bins : Dict[Tuple[str, Site], Dict[int, Set[str]]]
        each key is a tuple of a top branch and a site were it will be mutated,
        the associated value it's the set of bins.
    branch_bins_indices : Dict[Tuple[str, Site], BinPool]
        as ``branch_bins``, but the values are the indices of the bins that are
        available for the next mutation. The BinPool will refill each time
        it's emptied.
    epoch_bins_indices: BinPool
        indices of the bins that are available for the next mutation for any
        epoch and any site. The BinPool will refill each time it's emptied.
    """

    # Branch names + site are the keys in these 2 dicts.
    # We use them to better decide the next mutation on each branch.
    branch_bins: Dict[Tuple[str, Site], Dict[int, Set[str]]]
    branch_bins_indices: Dict[Tuple[str, Site], BinPool]
    epoch_bins_indices: BinPool

    def __init__(
        self,
        top_branches: Dict[str, Branch],
        bins: Dict[int, Set[str]],
        sites: List[Site],
    ) -> None:
        self.nbins = len(bins)
        self.epoch_bins_indices = BinPool.from_size(self.nbins)

        self.branch_bins = {}
        self.branch_bins_indices = {}
        for site in sites:
            for branch_name in top_branches.keys():
                self.branch_bins[(branch_name, site)] = deepcopy(bins)
                self.branch_bins_indices[(branch_name, site)] = BinPool.from_size(
                    self.nbins
                )

    def remove_from_memory(
        self, branch_name: str, site: Site, bin_idx: int, aa: str
    ) -> None:
        self.branch_bins_indices[(branch_name, site)].discard(bin_idx)
        self.epoch_bins_indices.discard(bin_idx)
        self.branch_bins[(branch_name, site)][bin_idx].discard(aa)

    def get_bin(
        self, branch_name: str, site: Site, old_aa: str, excluded_bins: Set[int]
    ) -> int:
        while True:
            bin_idx = self.__attempt_novel_bin_idx__(branch_name, site, excluded_bins)
            bin_aas = self.branch_bins[(branch_name, site)][bin_idx]

            bin_aas.discard(old_aa)
            if len(bin_aas) > 0:
                return bin_idx
            excluded_bins.add(bin_idx)
            if len(excluded_bins) == self.nbins:
                raise GenerationError(
                    f"Cannot mutate branch {branch_name}. "
                    "All amino acids from the available bins "
                    "have already been used."
                )

    def __attempt_novel_bin_idx__(
        self, branch_name: str, site: Site, excluded_bins: Set[int]
    ) -> int:
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
    """
    Used to generate a Dictionary with the top branches as keys and lists of
    mutations as values. It expects the list of Sites where the mutations
    will be carried out.

    Parameters
    ----------
    creation_config : Dict[str, Any]
        creation section from the input config file
    """

    use_bins: bool
    bins: Dict[int, Set[str]]
    aa_to_bin: Dict[str, Set[int]] = {"": set()}
    N_BINS: int
    N_SITES: int
    bins_criteria: str

    aa_distribution: Dict[str, float]

    memory: AminoAcidMemory

    def __init__(self, creation_config: Dict[str, Any]) -> None:
        self.bins = {}
        for i, each_bin in enumerate(creation_config["aa_bins"]):
            # Store bins
            self.bins[i] = set(each_bin)
            # Store in which bin is each amino acid
            for aa in each_bin:
                self.aa_to_bin[aa] = {i}
        self.N_BINS = len(self.bins)

        self.bins_criteria = creation_config["aa_bins_criteria"]
        self.aa_distribution = self.__initialize_probabilities__(creation_config)

    def __initialize_probabilities__(
        self, creation_config: Dict[str, Any]
    ) -> Dict[str, float]:
        if creation_config["aa_probability"] == "ReisBarletta":
            pd_reis_barletta_cdr = pd.read_csv(
                Path(Path(__file__).resolve().parent, "reis_barletta_cdr.csv")
            )
            aa_distribution = {
                r[1].AminoAcid: r[1].Probability
                for r in pd_reis_barletta_cdr.iterrows()
            }
        elif creation_config["aa_probability"] == "custom":
            aa_distribution = creation_config["custom"]
        else:
            pd_uniform = pd.read_csv(
                Path(Path(__file__).resolve().parent, "uniform.csv")
            )
            aa_distribution = {
                r[1].AminoAcid: r[1].Probability for r in pd_uniform.iterrows()
            }

        # Maybe I should do this check in cli.py
        bins_aas = set(chain.from_iterable(self.bins.values()))
        distro_aas = set(aa_distribution.keys())
        diff = set(bins_aas).symmetric_difference(distro_aas)
        if len(diff) != 0:
            raise GenerationError(
                "Fatal error. Amino acids from bins don't match the ones from the "
                f" probability distribution. Bins amino acids: {bins_aas}\n"
                f"Probability distribution amino acids: {distro_aas}\n"
                f"Different amino acids: {diff}."
            )
        return aa_distribution

    def __call__(
        self,
        top_branches: Dict[str, Branch],
        branches: int,
        sites: List[Site],
        logger: Logger,
    ) -> Dict[str, List[Mutation]]:
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
        # We'll use this set to make sure we're not creating identical branches.
        new_branch_names: Set[str] = set()

        while new_branches > 0:
            branch = top_branches[old_branches_pool.pop()]
            for site in sites:
                # No matter what, we'll decrease the branch counte,r so we're sure
                # it'll end. This means we may return less new branches than asked for.
                new_branches -= 1
                old_aa = branch.resnames[site.idx_chain][site.idx_residue]
                try:
                    new_aa = self.__select_aa__(branch.branch_name, site, old_aa)
                except GenerationError as e:
                    logger.warning(e)
                    continue

                mut = Mutation.from_site(site, old_aa=old_aa, new_aa=new_aa)
                branch_name, _ = branch.generate_name_resname(mut)
                if branch_name not in new_branch_names:
                    new_branch_names.add(branch_name)
                    # Add this mutation to this branch.
                    mutations[branch.branch_name].append(mut)

        return mutations

    def __select_aa__(self, branch_name: str, site: Site, old_aa: str) -> str:
        if self.bins_criteria == "without":
            excluded_bins: Set[int] = set(self.aa_to_bin[old_aa])
        else:
            # 'within'. All bins but the one that corresponds to the aa are excluded.
            excluded_bins: Set[int] = (
                set(range(0, len(self.bins))) - self.aa_to_bin[old_aa]
            )

        bin_idx = self.memory.get_bin(branch_name, site, old_aa, excluded_bins)
        # Normalize the probabilities of the AAs in the bin, so they add up to 1
        bin_aas = self.__get_norm_bin__(
            self.aa_distribution, self.memory.branch_bins[(branch_name, site)][bin_idx]
        )
        new_aa = choices(list(bin_aas.keys()), list(bin_aas.values()))[0]
        self.memory.remove_from_memory(branch_name, site, bin_idx, new_aa)

        return new_aa

    def __get_norm_bin__(
        self, distro: Dict[Union[str, int], float], aa_bin: Set[str]
    ) -> Dict[Union[str, int], float]:
        return self.__norm_distribution__(
            {aa: prob for aa, prob in distro.items() if aa in aa_bin}
        )

    @staticmethod
    def __norm_distribution__(
        distro: Dict[Union[str, int], float]
    ) -> Dict[Union[str, int], float]:
        tot = sum(distro.values())
        return {aa: prob / tot for aa, prob in distro.items()}
