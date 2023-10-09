from collections.abc import Iterable, ItemsView
from typing import List, Tuple, Dict, Set, Iterator, Final, Mapping, Any, Optional
from random import choice, sample
from collections import defaultdict
from logging import Logger

from attr import define, field

from interface import get_interfacing_residues
from locuaz.projectutils import Branch, Epoch
from locuaz.mutation import Mutation
from locuaz.abstractmutationgenerator import AbstractMutationGenerator


@define(frozen=True)
class Site:
    idx_chain: Final[int] = field(converter=int, kw_only=True)
    chainID: Final[str] = field(converter=str, kw_only=True)
    idx_residue: Final[int] = field(converter=int, kw_only=True)
    resSeq: Final[int] = field(converter=int, kw_only=True)


class SiteSelector:
    N_SITES: Final[int]
    excluded_pos: Set[int]
    probability: Final[str]
    use_interface: Final[bool]
    interface_probe_radius: Final[Optional[float]]

    idx_chain: int
    idx_residue: int

    def __init__(self, creation_config: Dict[str, Any], *, excluded_pos: Set[int]) -> None:
        self.N_SITES = creation_config["sites"]
        self.use_interface = creation_config["sites_interfacing"]
        self.interface_probe_radius = creation_config["sites_interfacing_probe_radius"]
        self.N_SITES = creation_config["sites"]
        self.excluded_pos = excluded_pos

    def __call__(self, epoch: Epoch, use_tleap: bool = False) -> List[Site]:
        # Get a branch to read the chainIDs and the resSeqs.
        try:
            any_branch = next(iter(epoch.top_branches.values()))
        except Exception:
            raise RuntimeError(f"No available branches on Epoch {epoch.id}. "
                               "It's likely that all of them failed during MD.")

        interface_resSeq: Set[int] = set()
        if self.use_interface:
            # Get the resSeq of all the residues in the interface and extend it
            interface_resSeq = get_interfacing_residues(any_branch.complex.pdb,
                                                        any_branch.chainIDs,
                                                        self.probe_radius,
                                                        use_tleap)
        else:
            # Pretend as if everything is on the interface.
            interface_resSeq = {resSeq for cdr in any_branch.resSeqs for resSeq in cdr}

        # Now, filter the mutating resSeqs.
        candidates_resSeq = [resSeq for cdr in any_branch.resSeqs for resSeq in cdr if
                             (resSeq not in self.excluded_pos) and (resSeq in interface_resSeq)]
        if len(candidates_resSeq) == 0:
            raise RuntimeError(f"No CDR residue on the interface that isn't excluded.\n"
                               f"Excluded positions: {self.excluded_pos}\n"
                               f"Interfacing positions: {self.interface_resSeq}.")
        sites: List[Site]
        for i in range(self.N_SITES):
            # Choose the position to mutate.
            mut_resSeq = choice(candidates_resSeq)
            # Now, get the remaining details associated to `mut_resSeq`, including the chain from where it came from
            for j, resSeqs in enumerate(any_branch.resSeqs):
                if mut_resSeq in resSeqs:
                    sites.append(Site(idx_chain=j,
                                      chainID=any_branch.chainIDs[j],
                                      idx_residue=resSeqs.index(mut_resSeq),
                                      resSeq=mut_resSeq))
        return sites


class AminoAcidSelector:
    use_bins: Final[bool]
    bins: Dict[str, List[str]] = {}
    N_BINS: Optional[Final[int]] = len(bins)
    bins_criteria: Final[str]
    aa_probability: Final[str]
    N_SITES: Final[int]

    remaining_bins: Set[int]
    excluded_aas: Set[str]

    def __init__(self, creation_config: Dict[str, Any], *, excluded_aas: Set[str],
                 use_tleap: bool = False) -> None:
        self.use_bins = creation_config["aa_bins_set"]
        self.bins = creation_config.get("bins", {})
        self.N_BINS = len(self.bins)
        self.remaining_bins = set(range(0, self.N_BINS))

        self.excluded_aas = excluded_aas

    def __call__(self, epoch: Epoch, branches: int, sites: Site) -> Dict[str, Mutation]:
        mutations = defaultdict(list)
        new_branches = branches
        remaining_branches = set(epoch.top_branches.keys())

        while new_branches != 0:
            for site in sites:
                branch = epoch.top_branches[remaining_branches.pop()]
                old_aa, new_aa = self.__get_random_aa__(branch, site.idx_chain, site.idx_residue)
                # Build the mutation object for this branch.
                mutation = Mutation(
                    chainID=site.chainID,
                    resSeq=site.resSeq,
                    old_aa=old_aa,
                    new_aa=new_aa,
                    chainID_idx=site.idx_chain,
                    resSeq_idx=site.idx_residue,
                )

                self.mutations[branch.branch_name].append(mutation)
                if len(remaining_branches) == 0:
                    # If all branches have been mutated at least once, and we still
                    # have branches to generate, restart `remaining_branches`.
                    remaining_branches = set(epoch.top_branches.keys())
                new_branches -= 1
        return mutations

    def __pop_random_aa__(self) -> str:
        cat_idx = choice(tuple(self.remaining_bins))
        self.remaining_bins.difference_update({cat_idx})

        shuffled_aas = sample(self.CAT_AAS[cat_idx], len(self.CAT_AAS[cat_idx]))
        for new_aa in shuffled_aas:
            if new_aa not in self.excluded_aas:
                self.excluded_aas.add(new_aa)
                return new_aa
        raise RuntimeError("Can't generate novel AA. This shouldn't happen.")

    def __get_random_aa__(self, branch: Branch, idx_chain: int, idx_residue: int) -> Tuple[str, str]:
        old_aa = branch.resnames[idx_chain][idx_residue]
        self.excluded_aas.add(old_aa)
        new_aa = self.__pop_random_aa__()

        if len(self.remaining_bins) == 0:
            # All categories have already been chosen from `N` times.
            # Allow all of them again for the `N+1` branch, except those
            # that are exhausted already.
            self.remaining_bins = set(range(self.N_BINS))
            for i in range(self.N_BINS):
                if set(self.CAT_AAS[i]).issubset(self.excluded_aas):
                    # All AAs from this category have already been chosen
                    self.remaining_bins.difference_update({i})
        return old_aa, new_aa


class MutationCreator(Mapping):
    """
    MutationCreator creates mutations by
        * choosing 1 or more positions after discarding the excluded ones and
        optionally discarding those that are not interfacing with the target.
        * optionally splits amino acids into user-defined bins. Then chooses the
        bin of the current residue ('within') or any of the other ones ('exclude').
        * within the potential amino acids, it chooses one of them following a
        uniform distribution ('uniform'), Reis&Barletta distribution ('ReisBarletta')
        or a custom user-defined probability ('custom')

        For each position, the chosen amino acid will be discarded from future
        consideration. This also applies to bins. If all bins were discarded and
        more mutations were asked for, the bins are reset but not the potential
        amino acids.
    """
    mutations: Dict[str, List[Mutation]]
    site_selector: SiteSelector
    aa_selector: AminoAcidSelector

    def __init__(
            self,
            epoch: Epoch,
            branches: int,
            creation_config: Dict[str, Any],
            *,
            excluded_aas: Set[str],
            excluded_pos: Set[int],
            use_tleap: bool = False
    ) -> None:
        self.site_selector = SiteSelector(creation_config, excluded_pos=excluded_pos, use_tleap=use_tleap)
        self.aa_selector = AminoAcidSelector(creation_config, excluded_aas=excluded_aas)

        sites = self.site_selector(epoch, use_tleap)
        self.mutations = self.aa_selector(epoch, branches, sites)

    def __getitem__(self, key: str) -> Iterable[Mutation]:
        return self.mutations[key]

    def __iter__(self) -> Iterator:
        return self.mutations.__iter__()

    def __contains__(self, value: Branch) -> bool:
        return self.mutations.__contains__(value)

    def __len__(self) -> int:
        return len(self.mutations)

    def __items__(self) -> ItemsView:
        return self.mutations.items()

    def __str__(self) -> str:
        return "Single Point Mutation from 4 amino acid groups"
