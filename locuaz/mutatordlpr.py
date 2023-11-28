import sys
from pathlib import Path
from typing import Optional, Union, Set, Tuple
import numpy as np

import MDAnalysis as mda
from Bio.SeqUtils import seq3

from locuaz.DLPacker.dlpacker import DLPacker
from locuaz.molecules import PDBStructure
from locuaz.mutation import Mutation
from locuaz.mutatordlp import MutatorDLPacker


class MutatorDLPackerReconstruct(MutatorDLPacker):
    radius: float

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
        super().__init__(bin_dir, radius, allowed_nonstandard_residues)
        self.radius = radius

    @staticmethod
    def __get_bridged_cys__(
        input_pdb: Union[PDBStructure, Path]
    ) -> Set[Tuple[int, str, str]]:
        u = mda.Universe(str(input_pdb))

        from MDAnalysis.analysis import distances

        cys_atoms = u.select_atoms("resname CYS and (name S or name SG)")
        dis_array = distances.self_distance_array(cys_atoms.positions)
        # Reshape the distances into a matrix
        n = len(cys_atoms)
        distancias = np.zeros((n, n))
        triu_indices = np.triu_indices_from(distancias, k=1)
        distancias[triu_indices] = dis_array
        distancias.T[triu_indices] = dis_array
        #
        for i in range(n):
            distancias[i, i] = sys.maxsize

        bridged_cys = set()
        for i, row in enumerate(distancias):
            j = np.argmin(row)
            if distancias[i, j] < 3.0:
                res_i = cys_atoms[i].residue
                bridged_cys.add((res_i.resnum, res_i.segment.segid, res_i.resname))
        return bridged_cys

    def __run__(
        self, input_pdb: Union[PDBStructure, Path], mutation: Mutation
    ) -> Tuple[PDBStructure, Set[int]]:
        input_pdb_fn = Path(input_pdb)

        dlp = DLPacker(
            str(input_pdb_fn),
            weights_path=Path(self.weights_path),
            lib_path=Path(self.lib_path),
            charges_path=Path(self.charges_path),
        )
        dlp.mutate_sequence(
            (mutation.resSeq, mutation.chainID, seq3(mutation.old_aa).upper()),
            seq3(mutation.new_aa).upper(),
        )

        # Get the surrounding residues
        targets = set(
            dlp.get_targets(
                target=(
                    mutation.resSeq,
                    mutation.chainID,
                    seq3(mutation.new_aa).upper(),
                ),
                radius=self.radius,
            )
        )
        # Don't move any cysteines that are forming bridges
        targets.difference_update(self.__get_bridged_cys__(input_pdb_fn))
        out_path = Path(input_pdb_fn.parent, "init_mutated.pdb")
        dlp.reconstruct_region(
            targets=list(targets), order="sequence", output_filename=str(out_path)
        )

        try:
            out_pdb = PDBStructure.from_path(out_path)
        except (FileNotFoundError, Exception) as e:
            raise f"{self} failed, not output mutated PDB. DLPacker object: {dlp}." from e

        return out_pdb, targets

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
            init_mutated_pdb, targets = self.__run__(nonwat_fix_pdb, mutation)
            excluded_res = " or ".join(
                f"(segid {chainID} and resnum {resSeq})"
                for resSeq, chainID, resname in targets
            )
            init_mutated_pdb = super().__fit_pdb__(
                init_mutated_pdb, nonwat_fix_pdb, excluded_res
            )
        except AssertionError as e:
            raise e
        dry_mut_pdb_fn = super().port_mutation(
            mutated_pdb=init_mutated_pdb,
            original_pdb=nonwat_pdb,
            mut=mutation,
            surrounding_residues=targets,
        )
        # Rejoin the mutated complex with water and ions
        overlapped_pdb_fn = Path(local_dir, "init_overlapped.pdb")
        overlapped_pdb = super().add_water(
            dry_mut_pdb_fn, wation_pdb, overlapped_pdb_fn
        )
        overlapped_pdb.set_cryst1_record(cryst1_record)

        return overlapped_pdb

    def __str__(self):
        return "MutatorDLPackerReconstruct"
