from attrs import define, field
from typing import List, Tuple

from projectutils import Iteration


@define(frozen=True)
class Mutation:
    chainID: str = field(converter=str, kw_only=True)
    resSeq: int = field(converter=int, kw_only=True)
    old_aa: str = field(converter=str, kw_only=True)
    new_aa: str = field(converter=str, kw_only=True)
    chainID_idx: int = field(converter=int, kw_only=True)
    resSeq_idx: int = field(converter=int, kw_only=True)

    # This is a horrible way to do this, I'll do it properly later.
    def new_name_resname(self, iteration: Iteration) -> Tuple[str, List[str]]:
        iter_name = ""
        new_iteration_resnames = []
        for idx, (chainID, resname) in enumerate(zip(iteration.chainIDs, iteration.resnames)):
            if idx == self.chainID_idx:
                # This is the mutated chainID
                new_resname = (
                        resname[: self.resSeq_idx]
                        + self.new_aa
                        + resname[self.resSeq_idx + 1:]
                )

            else:
                # This one remains the same
                new_resname = "".join([residue for residue in resname])
            iter_name += f"-{chainID}_{new_resname}"
            new_iteration_resnames.append(new_resname)

        # Drop the leading '-'
        return iter_name[1:], new_iteration_resnames

    def get_mda_sel(self, backbone_only=True) -> str:
        return f"segid {self.chainID} and resnum {self.resSeq}" + (" and backbone" if backbone_only else "")
