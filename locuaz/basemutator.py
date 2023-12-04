from abc import abstractmethod
from pathlib import Path
from typing import List, Iterable, Tuple, Optional, Union, Set
import warnings

import MDAnalysis as mda
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO

from locuaz.fileutils import FileHandle
from locuaz.projectutils import WorkProject, Epoch
from locuaz.molecules import PDBStructure
from locuaz.mutation import Mutation
from locuaz.primitives import AA_MAP


class BaseMutator:
    """BaseMutator offers 4 static methods to help development of new mutators: split_solute_solvent(), fix_pdb(),
    port_mutation() and add_water().
    A Mutator may output a PDB that has waters overlapping with the newly mutated residue.
    Parameters
    ----------
        bin_dir : Path
            path to the binary file.
        radius : float
            radius freesasa will use to define the Ag-Ab interface. Default: 5A.

    """

    @abstractmethod
    def __init__(
        self,
        bin_dir: Path,
        radius: float = 5,
        allowed_nonstandard_residues: Optional[Set[str]] = None,
    ) -> None:
        pass

    def on_pdb(
        self,
        input_pdb: PDBStructure,
        local_dir: Path,
        *,
        mutation: Mutation,
        selection_complex: Optional[str] = None,
        selection_wations: Optional[str] = None,
    ) -> PDBStructure:
        pass

    @staticmethod
    def split_solute_solvent(
        pdb_in: Union[PDBStructure, Path],
        *,
        selection_complex: Optional[str] = None,
        selection_wations: Optional[str] = None,
    ) -> Tuple[PDBStructure, PDBStructure]:
        pdb_in_path = Path(pdb_in)
        u = mda.Universe(str(pdb_in_path))

        # Protein
        nonwat_pdb_fn = Path(pdb_in_path.parent, "init_nonwat.pdb")
        if selection_complex:
            u.atoms.select_atoms(selection_complex).write(str(nonwat_pdb_fn))  # type: ignore
        else:
            u.atoms.select_atoms("protein").write(str(nonwat_pdb_fn))  # type: ignore

        # Water and ions
        wation_pdb_fn = Path(pdb_in_path.parent, "init_wation.pdb")
        if selection_wations:
            u.atoms.select_atoms(selection_wations).write(str(wation_pdb_fn))  # type: ignore
        else:
            u.atoms.select_atoms("not protein").write(str(wation_pdb_fn))  # type: ignore

        return PDBStructure.from_path(nonwat_pdb_fn), PDBStructure.from_path(
            wation_pdb_fn
        )

    def fix_pdb(self, pdb_in: Union[Path, PDBStructure]) -> Path:
        """
        Uses AA_MAP to map non-standard amino acids to standard ones (eg: 'CY2' to 'CYS')
        and then removes non-amino acidic molecules (eg: 'ZN' and other ligands).
        This is to prevent mutators from choking on some proteins.

        Parameters
        ----------
        pdb_in : Union[Path, PDBStructure]
            input dried PDB without water nor ions which may contain
            non-standard residues and non-amino acidic molecules

        Returns
        -------
        the output PDB : Path
            PDB without non-amino acidic molecules and with standardized resnames.
        """
        pdb_in_fn = Path(pdb_in)
        u = mda.Universe(str(pdb_in_fn))
        not_prot: List[str] = []
        for res in u.residues:  # type: ignore
            try:
                res.resname = AA_MAP[res.resname]
            except KeyError:
                if res.resname not in self.allowed_nonstandard_residues:
                    # Not even an amino acid, removing it.
                    not_prot.append(res.resname)
        not_prot_sel = " and ".join(f"not resname {res}" for res in not_prot)

        pdb_out_fn = Path(pdb_in_fn.parent, "init_nonwat_fix.pdb")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if not_prot_sel == "":
                # No non-protein residues
                u.atoms.write(str(pdb_out_fn))
            else:
                u.select_atoms(not_prot_sel).write(str(pdb_out_fn))

        return pdb_out_fn

    @staticmethod
    def add_water(
        solute_path: Union[Path, PDBStructure],
        solvent_path: Union[Path, PDBStructure],
        out_path: Path,
    ) -> PDBStructure:
        u = mda.Universe(str(solute_path))
        v = mda.Universe(str(solvent_path))

        # TODO: get a `charge_diff` parameter, if not 0, count NA/Na and CL/Cl ions from the solvent
        # TODO: 1st try to remove an ion to match the `charge_diff`, if impossible, then
        # TODO: add `charge_diff` NA/Na or CL/Cl ions by using the `ion_contour` code
        # TODO: (using mda.Merge to get proper numbering) this code should be in BaseMutator.

        # Prevent MDA from creating its own chainID like 'X' or 'SYSTEM'
        v.atoms.segments.segids = " "  # type: ignore
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mda.Merge(u.atoms, v.atoms).atoms.write(str(out_path))  # type: ignore

        return PDBStructure.from_path(out_path)

    @staticmethod
    def port_mutation(
        mutated_pdb: Union[Path, PDBStructure],
        original_pdb: Union[Path, PDBStructure],
        mut: Mutation,
        surrounding_residues: Optional[Set[Tuple[int, str, str]]] = None,
    ) -> Path:
        mutated_pdb_fn = Path(mutated_pdb)
        parsero = PDBParser(QUIET=True)
        mut_pdb = parsero.get_structure("foo", file=mutated_pdb_fn)
        orig_pdb = parsero.get_structure("bar", file=Path(original_pdb))

        new_model = mut_pdb[0]
        model = orig_pdb[0]

        if surrounding_residues:
            for resSeq, chainID, resname in surrounding_residues:
                new_segment = new_model[chainID]
                try:
                    new_aa = new_segment.child_dict[(" ", resSeq, " ")]
                except KeyError:
                    # if the residue is not an amino acid, BioPDB names it this way.
                    new_aa = new_segment.child_dict[(f"H_{resname}", resSeq, " ")]

                segment = model[chainID]
                first_resSeq = segment.child_list[0].id[1]
                segment.child_list[resSeq - first_resSeq] = new_aa
        else:
            new_segment = new_model[mut.chainID]
            new_aa = new_segment.child_dict[(" ", mut.resSeq, " ")]

            segment = model[mut.chainID]
            first_resSeq = segment.child_list[0].id[1]
            segment.child_list[mut.resSeq - first_resSeq] = new_aa

        io = PDBIO()
        io.set_structure(orig_pdb)
        out_path = Path(mutated_pdb_fn.parent, "init_nonwat_mut.pdb")
        io.save(str(out_path))

        return out_path

    def assert_outfile(
        self,
        out_file: Union[str, Path, FileHandle],
        *,
        stdout: str,
        stderr: str,
        command: str,
    ) -> Path:
        out_file_path = Path(out_file)
        assert out_file_path.is_file(), f"""{self} error. Can't parse: {out_file_path}
from:
{command}
with stdout:
{stdout}
and stderr:
{stderr}
"""
        return out_file_path


def memorize_mutations(
    work_pjct: WorkProject, new_epoch: Epoch, mutations: Iterable[Mutation]
) -> None:
    if work_pjct.has_memory:
        new_epoch.mutated_positions = {mutation.resSeq for mutation in mutations}
        work_pjct.mutated_positions.appendleft(new_epoch.mutated_positions)
    return
