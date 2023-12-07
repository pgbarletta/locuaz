import shutil as sh
import subprocess as sp
from pathlib import Path
from typing import Optional, Union

from locuaz.fileutils import FileHandle
from locuaz.molecules import PDBStructure
from locuaz.mutation import Mutation
from locuaz.basemutator import BaseMutator


class MutatorEvoEF2(BaseMutator):
    bin_path: FileHandle

    def __init__(self, bin_dir: Path, radius: float = 5) -> None:
        self.bin_path = FileHandle(Path(bin_dir, "evoef2"))

    @staticmethod
    def __evoef2_file__(mut: Mutation, local_dir: Path) -> FileHandle:
        """__evoef2_file__ creates a .txt file for EvoEF2's --command=BuildMutant

        Returns:
            FileHandle: .txt file with mutations as: \
                WT_AA_OneLetterCode Chain Resnum MUT_AA_OneLetterCode (no whitespace)
        """
        mut_string = f"{mut.old_aa}{mut.chainID}{mut.resSeq}{mut.new_aa};"
        mutlist_path = local_dir / f"mutlist_evoef2_{mut_string[:-1]}.txt"

        with open(mutlist_path, "w") as file:
            file.write(mut_string)

        return FileHandle(mutlist_path)

    def __run__(
            self, input_pdb: Union[PDBStructure, Path], mutation: Mutation
    ) -> PDBStructure:

        input_pdb_fn = Path(input_pdb)
        local_dir = input_pdb_fn.parent

        mutlist = self.__evoef2_file__(mutation, local_dir)

        # Using relative paths to keep them short.
        comando_evoef2 = f"{self.bin_path} --command=BuildMutant --pdb={input_pdb_fn.name} --mutant_file={mutlist.name_ext()}"

        p = sp.run(
            comando_evoef2,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            cwd=local_dir,
            shell=True,
            text=True,
        )
        # This is EvoEF2's naming convention
        mut_path = Path(local_dir, f"{input_pdb_fn.stem}_Model_0001.pdb")
        try:
            self.assert_outfile(mut_path, stdout=p.stdout, stderr=p.stderr, command=comando_evoef2)
        except AssertionError as e:
            raise e
        out_path = Path(local_dir, "init_mutated.pdb")

        sh.move(mut_path, out_path)

        return PDBStructure.from_path(out_path)

    def on_pdb(
            self,
            input_pdb: PDBStructure,
            local_dir: Path,
            *,
            mutation: Mutation,
            selection_complex: Optional[str] = None,
            selection_wations: Optional[str] = None,
    ) -> PDBStructure:
        # Get the system's box size after the NPT run, to add it later onto the
        # mutated PDB system. The PDB format has less precision for the box parameters
        # than the GRO format, so there may be a difference in the last digit for the
        # lengths (eg: 12.27215 to 12.27210) and the angles (6.13607 to 6.13605).
        cryst1_record = input_pdb.get_cryst1_record()
        wt_pdb_fn = Path(input_pdb)

        nonwat_pdb, wation_pdb = super().split_solute_solvent(
            wt_pdb_fn,
            selection_complex=selection_complex,
            selection_wations=selection_wations,
        )

        nonwat_fix_pdb = super().fix_pdb(nonwat_pdb)
        try:
            init_mutated_pdb = self.__run__(nonwat_fix_pdb, mutation)
        except AssertionError as e:
            raise e
        dry_mut_pdb_fn = super().port_mutation(mutated_pdb=init_mutated_pdb, original_pdb=nonwat_pdb, mut=mutation)
        # Rejoin the mutated complex with water and ions
        overlapped_pdb_fn = Path(local_dir, "init_overlapped.pdb")
        overlapped_pdb = super().add_water(
            dry_mut_pdb_fn, wation_pdb, overlapped_pdb_fn
        )
        overlapped_pdb.set_cryst1_record(cryst1_record)

        return overlapped_pdb

    def __str__(self):
        return "MutatorEvoEF2"
