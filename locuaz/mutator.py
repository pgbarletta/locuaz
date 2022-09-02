from pathlib import Path
from abc import ABCMeta, abstractmethod
from attrs import define, field
from typing import Dict, List, Tuple
import subprocess as sp
from Bio.SeqUtils import seq3
from projectutils import Iteration
from molecules import PDBStructure
from biobb_analysis.gromacs.gmx_trjconv_str import GMXTrjConvStr
from biobb_model.model.mutate import Mutate
from primitives import launch_biobb
from fileutils import FileHandle


@define(frozen=True)
class Mutation:

    chainID: str = field(converter=str, kw_only=True)
    resSeq: int = field(converter=int, kw_only=True)
    old_aa: str = field(converter=str, kw_only=True)
    new_aa: str = field(converter=str, kw_only=True)
    chainID_idx: int = field(converter=int, kw_only=True)
    resSeq_idx: int = field(converter=int, kw_only=True)

    def biobb_string(self) -> str:
        """biobb_string returns a string for biobb.model.mutate.mutate

        Returns:
            str: Chain:WT_AA_ThreeLeterCode Resnum MUT_AA_ThreeLeterCode (no whitespace)
        """
        return f"{self.chainID}:{seq3(self.old_aa)}{self.resSeq}{seq3(self.new_aa)}"

    def evoef2_file(self, output_path: Path) -> FileHandle:
        """evoef2_file returns a .txt file for EvoEF2's --command=BuildMutant

        Returns:
            FileHandle: .txt file with mutations as: \
                WT_AA_ThreeLeterCode Chain Resnum MUT_AA_ThreeLeterCode (no whitespace)
        """
        mut_string = f"{self.old_aa}{self.chainID}{self.resSeq}{self.new_aa};"
        with open(output_path, "w") as file:
            file.write(mut_string)

        return FileHandle(output_path)

    def new_name_resname(self, iter: Iteration) -> Tuple[str, List[str]]:
        iter_name = ""
        new_iteration_resnames = []
        for chainID, resname in zip(iter.chainIDs, iter.resnames):
            if chainID == self.chainID:
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
    def __init__(self):
        pass

    @abstractmethod
    def __call__(self, input_pdb: PDBStructure, mutation: Mutation) -> PDBStructure:
        pass


class MutatorEvoEF2(AbstractMutator):
    bin_path: FileHandle

    def __init__(self, config: Dict):
        self.bin_path = FileHandle(
            Path(config["paths"]["scoring_functions"]) / "evoef2/EvoEF2"
        )

    def __call__(self, input_pdb: PDBStructure, mutation: Mutation) -> PDBStructure:

        wrk_dir = input_pdb.file.path.parent
        input_mutlist_fn = wrk_dir / "mutlist.txt"
        mutation.evoef2_file(input_mutlist_fn)

        comando_evoef2 = (
            str(self.bin_path)
            + " --command=BuildMutant"
            + " --pdb="
            + str(input_pdb.file.path)
            + " --mutant_file="
            + str(input_mutlist_fn)
        )

        sp.run(
            comando_evoef2,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            cwd=wrk_dir,
            shell=True,
            text=True,
        )
        # This is EvoEF2's naming convention
        temp_out_pdb = wrk_dir / (
            input_pdb.file.path.name.split(".")[0] + "_Model_0001.pdb"
        )

        input_mutlist_fn.unlink()

        return PDBStructure.from_path(temp_out_pdb)


class MutatorBiobb(AbstractMutator):
    def __init__(self, config: Dict):
        self.bin_path = FileHandle(
            Path(config["paths"]["scoring_functions"]) / "evoef2/EvoEF2"
        )

    def __call__(self, input_pdb: PDBStructure, mutation: Mutation) -> PDBStructure:

        wrk_dir = input_pdb.file.path.parent
        temp_out_pdb = wrk_dir / (
            input_pdb.file.path.name.split(".")[0] + "_Model_0001.pdb"
        )
        biobb_mutador = Mutate(
            input_pdb_path=str(input_pdb.file.path),
            output_pdb_path=str(temp_out_pdb),
            properties={"mutation_list": mutation.biobb_string()},
        )
        launch_biobb(biobb_mutador)

        return PDBStructure.from_path(temp_out_pdb)
