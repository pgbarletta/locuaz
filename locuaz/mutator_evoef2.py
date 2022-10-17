from pathlib import Path
import subprocess as sp

from molecules import PDBStructure
from fileutils import FileHandle
from mutator import AbstractMutator, Mutation


class MutatorEvoEF2(AbstractMutator):
    bin_path: FileHandle

    def __init__(self, bin_dir: Path) -> None:
        self.bin_path = FileHandle(Path(bin_dir, "EvoEF2"))

    def __evoef2_file__(self, mut: Mutation, output_path: Path) -> FileHandle:
        """__evoef2_file__ creates a .txt file for EvoEF2's --command=BuildMutant

        Returns:
            FileHandle: .txt file with mutations as: \
                WT_AA_ThreeLeterCode Chain Resnum MUT_AA_ThreeLeterCode (no whitespace)
        """
        mut_string = f"{mut.old_aa}{mut.chainID}{mut.resSeq}{mut.new_aa};"
        with open(output_path, "w") as file:
            file.write(mut_string)

        return FileHandle(output_path)

    def __call__(self, input_pdb: PDBStructure, mutation: Mutation) -> PDBStructure:

        wrk_dir = input_pdb.file.path.parent
        input_mutlist_fn = wrk_dir / "mutlist.txt"
        self.__evoef2_file__(mutation, input_mutlist_fn)

        # Using relative paths to keep them short.
        comando_evoef2 = (
            f"{self.bin_path} --command=BuildMutant --pdb={input_pdb.file.path.name} "
            "--mutant_file=mutlist.txt"
            # + str(input_mutlist_fn.name)
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
