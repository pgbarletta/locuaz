from pathlib import Path
from typing import Optional, Union, Set

import MDAnalysis as mda
from MDAnalysis.analysis import align
from Bio.SeqUtils import seq3

from locuaz.DLPacker.dlpacker import DLPacker
from locuaz.fileutils import FileHandle
from locuaz.molecules import PDBStructure
from locuaz.mutation import Mutation
from locuaz.basemutator import BaseMutator


class MutatorDLPacker(BaseMutator):
    weights_path: FileHandle
    lib_path: FileHandle
    charges_path: FileHandle
    allowed_nonstandard_residues: Set[str]

    def __init__(
        self,
        bin_dir: Path,
        radius: float = 5,
        allowed_nonstandard_residues: Optional[Set[str]] = None,
    ) -> None:
        """

        Parameters
        ----------
        bin_dir
        radius
        allowed_nonstandard_residues : Optional[set]
                Any residue not present in AA_MAP or in this set will be discarded
                from the output PDB. Useful when mutating complexes with ligands as
                antigens.
        """
        self.weights_path = FileHandle(Path(bin_dir, "DLPacker_weights.h5"))
        self.lib_path = FileHandle(Path(bin_dir, "library.npz"))
        self.charges_path = FileHandle(Path(bin_dir, "charges.rtp"))
        if allowed_nonstandard_residues:
            self.allowed_nonstandard_residues = allowed_nonstandard_residues
        else:
            self.allowed_nonstandard_residues = set()

    def __fit_pdb__(
        self,
        mutated_pdb: Union[Path, PDBStructure],
        wt_pdb: Union[Path, PDBStructure],
        excluded_selection: str,
    ) -> Path:
        u = mda.Universe(str(mutated_pdb))
        v = mda.Universe(str(wt_pdb))

        align.alignto(u.atoms, v.atoms, select="backbone")
        align.alignto(
            u.atoms, v.atoms, select=f"backbone and not ({excluded_selection})"
        )

        mutated_pdb_fn = Path(mutated_pdb)
        out_pdb = Path(
            mutated_pdb_fn.parent, mutated_pdb_fn.stem + "_fit" + mutated_pdb_fn.suffix
        )
        u.atoms.write(str(out_pdb))

        return out_pdb

    def __run__(
        self, input_pdb: Union[PDBStructure, Path], mutation: Mutation
    ) -> PDBStructure:
        input_pdb_fn = Path(input_pdb)

        dlp = DLPacker(
            str(input_pdb_fn),
            weights_path=Path(self.weights_path),
            lib_path=Path(self.lib_path),
            charges_path=Path(self.charges_path),
        )
        dlp.mutate_residue(
            (mutation.resSeq, mutation.chainID, seq3(mutation.old_aa).upper()),
            seq3(mutation.new_aa).upper(),
        )

        out_path = Path(input_pdb_fn.parent, "init_mutated.pdb")
        dlp.save_structure(str(out_path))

        try:
            out_pdb = PDBStructure.from_path(out_path)
        except (FileNotFoundError, Exception) as e:
            raise f"{self} failed, not output mutated PDB. DLPacker object: {dlp}." from e

        return out_pdb

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
            init_mutated_pdb = self.__fit_pdb__(
                init_mutated_pdb, nonwat_fix_pdb, mutation.get_mda_sel()
            )
        except AssertionError as e:
            raise e
        dry_mut_pdb_fn = super().port_mutation(
            mutated_pdb=init_mutated_pdb, original_pdb=nonwat_pdb, mut=mutation
        )
        # Rejoin the mutated complex with water and ions
        overlapped_pdb_fn = Path(local_dir, "init_overlapped.pdb")
        overlapped_pdb = super().add_water(
            dry_mut_pdb_fn, wation_pdb, overlapped_pdb_fn
        )
        overlapped_pdb.set_cryst1_record(cryst1_record)

        return overlapped_pdb

    def __str__(self):
        return "MutatorDLPacker"
