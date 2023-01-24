from pathlib import Path
from abc import ABCMeta, abstractmethod
from attrs import define, field
from typing import Iterable, List, Iterable, Tuple, Optional


from projectutils import Iteration, WorkProject, Epoch
from molecules import PDBStructure
from complex import GROComplex


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
    def __run__(self, input_pdb: PDBStructure, mutation: Mutation) -> PDBStructure:
        pass

    @abstractmethod
    def on_pdb(
        self,
        input_pdb: PDBStructure,
        local_dir: Path,
        *,
        mutation: Mutation,
        selection_protein: Optional[str] = None,
        selection_wations: Optional[str] = None,
    ) -> PDBStructure:
        pass

    # @abstractmethod
    # def on_grocomplex(
    #     self, cpx: GROComplex, *, mutation: Mutation, split_solvent=True
    # ) -> PDBStructure:
    #     pass


def memorize_mutations(
    work_pjct: WorkProject, new_epoch: Epoch, mutations: Iterable[Mutation]
) -> None:

    mutated_aminoacids = [mutation.new_aa for mutation in mutations]
    mutated_positions = [mutation.resSeq for mutation in mutations]
    new_epoch.mutated_positions = set(mutated_positions)

    if not work_pjct.has_memory:
        return
    else:
        # TODO: memorize amino acid as well.
        work_pjct.mutated_positions.appendleft(new_epoch.mutated_positions)
