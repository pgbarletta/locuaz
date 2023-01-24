from pathlib import Path
import subprocess as sp
from typing import Optional, Union, List
import shutil as sh

import MDAnalysis as mda
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO

from molecules import PDBStructure
from primitives import AA_MAP
from fileutils import FileHandle, catenate_pdbs
from mutator import AbstractMutator, Mutation
from molutils import split_solute_solvent
from amberutils import fix_pdb


class MutatorEvoEF2(AbstractMutator):
    bin_path: FileHandle

    def __init__(self, bin_dir: Path) -> None:
        self.bin_path = FileHandle(Path(bin_dir, "EvoEF2"))

    def __evoef2_file__(self, mut: Mutation, local_dir: Path) -> FileHandle:
        """__evoef2_file__ creates a .txt file for EvoEF2's --command=BuildMutant

        Returns:
            FileHandle: .txt file with mutations as: \
                WT_AA_OneLetterCode Chain Resnum MUT_AA_OneLetterCode (no whitespace)
        """
        mut_string = f"{mut.old_aa}{mut.chainID}{mut.resSeq}{mut.new_aa};"
        mutlist_path = local_dir / f"mutlist_{mut_string[:-1]}.txt"

        with open(mutlist_path, "w") as file:
            file.write(mut_string)

        return FileHandle(mutlist_path)

    def __run__(
        self, input_pdb: Union[PDBStructure, Path], mutation: Mutation
    ) -> PDBStructure:

        input_pdb_fn = Path(input_pdb)
        local_dir = input_pdb_fn.parent

        mutlist = self.__evoef2_file__(mutation, local_dir)

        # premut = Path(local_dir, "premut.pdb")
        # fix_gromacs_pdb(input_pdb.file.path, premut)

        # Using relative paths to keep them short.
        comando_evoef2 = f"{self.bin_path} --command=BuildMutant --pdb={input_pdb_fn.name} --mutant_file={mutlist}"

        sp.run(
            comando_evoef2,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            cwd=local_dir,
            shell=True,
            text=True,
        )
        # This is EvoEF2's naming convention
        mut_path = Path(local_dir, f"{input_pdb_fn.stem}_Model_0001.pdb")
        out_path = Path(local_dir, "init_mutated.pdb")

        sh.move(mut_path, out_path)

        return PDBStructure.from_path(out_path)

    def on_pdb(
        self,
        input_pdb: PDBStructure,
        local_dir: Path,
        *,
        mutation: Mutation,
        selection_protein: Optional[str] = None,
        selection_wations: Optional[str] = None,
    ) -> PDBStructure:
        # Get the system's box size after the NPT run, to add it later onto the
        # mutated PDB system. The PDB format has less precision for the box parameters
        # than the GRO format, so there may be a difference in the last digit for the
        # lengths (eg: 12.27215 to 12.27210) and the angles (6.13607 to 6.13605).
        # That's why GROComplex.from_pdb() also uses editconf.
        cryst1_record = input_pdb.get_cryst1_record()

        # Copy starting PDB into the new iteration folder (local_dir).
        wt_pdb_fn = Path(local_dir, "init_wt.pdb")
        sh.copy(Path(input_pdb), wt_pdb_fn)

        nonwat_pdb, wation_pdb = split_solute_solvent(
            wt_pdb_fn,
            selection_protein=selection_protein,
            selection_wations=selection_wations,
        )

        nonwat_fix_pdb = self.__fix_pdb__(nonwat_pdb)
        init_mutated_pdb = self.__run__(nonwat_fix_pdb, mutation)

        dry_mut_pdb_fn = self.__port_mutation__(
            mutated_pdb=init_mutated_pdb,
            original_pdb=nonwat_pdb,
            mut=mutation,
        )

        # Rejoin the mutated complex with water and ions
        overlapped_pdb_fn = Path(local_dir, "init_overlapped.pdb")
        overlapped_pdb = self.__add_water__(
            dry_mut_pdb_fn, wation_pdb, overlapped_pdb_fn
        )
        overlapped_pdb.set_cryst1_record(cryst1_record)

        return overlapped_pdb

    def __port_mutation__(
        self,
        *,
        mutated_pdb: Union[Path, PDBStructure],
        original_pdb: Union[Path, PDBStructure],
        mut: Mutation,
    ) -> Path:
        mutated_pdb_fn = Path(mutated_pdb)
        parsero = PDBParser(QUIET=True)
        mut_pdb = parsero.get_structure("foo", file=mutated_pdb_fn)
        orig_pdb = parsero.get_structure("bar", file=Path(original_pdb))

        new_model = mut_pdb[0]
        new_segment = new_model[mut.chainID]
        new_aa = new_segment.child_dict[(" ", mut.resSeq, " ")]

        model = orig_pdb[0]
        segment = model[mut.chainID]
        first_resSeq = segment.child_list[0].id[1]
        segment.child_list[mut.resSeq - first_resSeq] = new_aa

        io = PDBIO()
        io.set_structure(orig_pdb)
        out_path = Path(mutated_pdb_fn.parent, "init_nonwat_mut.pdb")
        io.save(str(out_path))

        return out_path

    def __add_water__(
        self,
        solute_path: Union[Path, PDBStructure],
        solvent_path: Union[Path, PDBStructure],
        out_path: Path,
    ) -> PDBStructure:
        u = mda.Universe(str(solute_path))
        v = mda.Universe(str(solvent_path))

        # Prevent MDA from creating its own chainID like 'X' or 'SYSTEM'
        v.atoms.segments.segids = " "  # type: ignore
        mda.Merge(u.atoms, v.atoms).atoms.write(str(out_path))  # type: ignore

        return PDBStructure.from_path(out_path)

    def __fix_pdb__(self, pdb_in: Union[Path, PDBStructure]) -> Path:
        pdb_in_fn = Path(pdb_in)
        u = mda.Universe(str(pdb_in_fn))
        not_prot: List[str] = []
        for res in u.residues:  # type: ignore
            try:
                res.resname = AA_MAP[res.resname]
            except KeyError:
                not_prot.append(res.resname)
        sel = "not resname " + " not resname and ".join([res for res in not_prot])

        pdb_out_fn = Path(pdb_in_fn.parent, "init_nonwat_fix.pdb")
        u.select_atoms(sel).write(str(pdb_out_fn))

        return pdb_out_fn
