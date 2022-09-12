from pathlib import Path
from abc import ABCMeta, abstractmethod
from attrs import define, field
from typing import List, Tuple


from projectutils import Iteration
from molecules import PDBStructure
from fileutils import FileHandle


@define(frozen=True)
class Mutation:

    chainID: str = field(converter=str, kw_only=True)
    resSeq: int = field(converter=int, kw_only=True)
    old_aa: str = field(converter=str, kw_only=True)
    new_aa: str = field(converter=str, kw_only=True)
    chainID_idx: int = field(converter=int, kw_only=True)
    resSeq_idx: int = field(converter=int, kw_only=True)

    # This is a horrible way to do this, I'll do it properly later.
    def new_name_resname(self, iter: Iteration) -> Tuple[str, List[str]]:
        iter_name = ""
        new_iteration_resnames = []
        for idx, (chainID, resname) in enumerate(zip(iter.chainIDs, iter.resnames)):
            if idx == self.chainID_idx:
                # This is the mutated chainID
                new_resname = (
                    resname[: self.resSeq_idx]
                    + self.new_aa
                    + resname[self.resSeq_idx + 1 :]
                )

            else:
                # This one remains the same
                new_resname = "".join([residue for residue in resname])
            iter_name += f"-{chainID}_{new_resname}"
            new_iteration_resnames.append(new_resname)

        # Drop the leading '-'
        return iter_name[1:], new_iteration_resnames


class AbstractMutator(metaclass=ABCMeta):
    @abstractmethod
    def __init__(self, bin_dir: Path) -> None:
        pass

    @abstractmethod
    def __call__(self, input_pdb: PDBStructure, mutation: Mutation) -> PDBStructure:
        pass
