from collections.abc import Iterable, ItemsView
from logging import Logger
from typing import List, Dict, Set, Iterator, Mapping, Any, Optional

from locuaz.aminoacidselector import AminoAcidSelector
from locuaz.mutation import Mutation
from locuaz.projectutils import Branch, Epoch
from locuaz.siteselector import SiteSelector


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
            top_branches: Dict[str, Branch],
            branches: int,
            creation_config: Dict[str, Any],
            *,
            excluded_sites: Set[int],
            amber_numbering: bool = False,
            logger: Logger
    ) -> None:
        any_branch = next(iter(top_branches.values()))

        self.site_selector = SiteSelector(any_branch.resSeqs,
                                          any_branch.chainIDs,
                                          excluded_sites,
                                          creation_config,
                                          amber_numbering=amber_numbering)
        self.aa_selector = AminoAcidSelector(creation_config)

        sites = self.site_selector(top_branches, logger)
        self.mutations = self.aa_selector(top_branches, branches, sites, logger)

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
