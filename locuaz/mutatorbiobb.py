from pathlib import Path

from Bio.SeqUtils import seq3
from biobb_model.model.mutate import Mutate

from locuaz.molecules import PDBStructure
from locuaz.primitives import launch_biobb
from locuaz.basemutator import BaseMutator
from locuaz.mutation import Mutation


class MutatorBiobb(BaseMutator):
    def __init__(self, bin_dir: Path) -> None:
        pass

    def __biobb_string__(self, mut: Mutation) -> str:
        """biobb_string returns a string for biobb.model.mutate.mutate

        Returns:
            str: Chain:WT_AA_ThreeLeterCode Resnum MUT_AA_ThreeLeterCode (no whitespace)
        """
        return f"{mut.chainID}:{seq3(mut.old_aa)}{mut.resSeq}{seq3(mut.new_aa)}"

    def __call__(self, input_pdb: PDBStructure, mutation: Mutation) -> PDBStructure:

        wrk_dir = input_pdb.file.path.parent
        # Using EvoEF2's naming convention
        temp_out_pdb = wrk_dir / (
            input_pdb.file.path.name.split(".")[0] + "_Model_0001.pdb"
        )
        biobb_mutador = Mutate(
            input_pdb_path=str(input_pdb.file.path),
            output_pdb_path=str(temp_out_pdb),
            properties={"mutation_list": self.__biobb_string__(mutation)},
        )
        launch_biobb(biobb_mutador)

        return PDBStructure.from_path(temp_out_pdb)
