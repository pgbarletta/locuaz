import shutil as sh
import sys
import warnings
from logging import Logger
from pathlib import Path

from locuaz.amberutils import fix_pdb
from locuaz.complex import GROComplex
from locuaz.fileutils import DirHandle
from locuaz.gromacsutils import remove_overlapping_solvent
from locuaz.molecules import PDBStructure
from locuaz.mutationgenerators import mutation_generators
from locuaz.basemutator import memorize_mutations
from locuaz.mutators import mutators
from locuaz.projectutils import WorkProject, Epoch, Branch


def initialize_new_epoch(work_pjct: WorkProject, log: Logger) -> None:
    """initialize_new_epoch(): This is a specific protocol, others will be added

    Args:
        work_pjct (WorkProject):
        log (Logger):
    """
    name = work_pjct.config["main"]["name"]
    old_epoch = work_pjct.epochs[-1]
    epoch_id = old_epoch.id + 1
    current_epoch = Epoch(epoch_id, branches={}, nvt_done=False, npt_done=False)

    # Create required mutator
    mutator = mutators[work_pjct.config["mutation"]["mutator"]](
        work_pjct.config["paths"]["mutator"], work_pjct.config["mutation"]["reconstruct_radius"]
    )
    # Create required mutation generator and generate mutation.
    generator = mutation_generators[work_pjct.config["generation"]["generator"]]
    successful_mutations = 0

    if work_pjct.config["protocol"]["constant_width"]:
        branches = work_pjct.config["protocol"]["new_branches"]
    else:
        branches = work_pjct.config["protocol"]["new_branches"] * len(old_epoch.top_branches)

    # Usually, this `while` would only be executed once, unless the Mutator program fails to perform a mutation.
    while successful_mutations < branches:
        mutation_generator = generator(
            old_epoch,
            branches - successful_mutations,
            excluded_aas=work_pjct.get_mem_aminoacids(),
            excluded_pos=work_pjct.get_mem_positions(),
            use_tleap=work_pjct.config["md"]["use_tleap"],
            logger=log,
            probe_radius=work_pjct.config["generation"]["probe_radius"])

        for old_branch_name, mutations in mutation_generator.items():
            old_iter = old_epoch.top_branches[old_branch_name]

            if work_pjct.config["md"]["use_tleap"]:
                # GROMACS renumbers resSeqs to strided numbering. If using Amber's continuous
                # numbering, this will result in the wrong mutating_resSeq.
                # Backup the PDB before running pdb4amber
                pdb_path = Path(old_iter.complex.pdb)
                pre_fix_pdb = Path(old_iter.dir_handle, f"preAmberPDBFixer_{pdb_path.stem}.pdb")
                sh.move(pdb_path, pre_fix_pdb)

                old_pdb = fix_pdb(pre_fix_pdb, pdb_path)
            else:
                old_pdb = old_iter.complex.pdb

            for mutation in mutations:
                branch_name, branch_resnames = old_iter.generate_name_resname(mutation)
                branch_path = Path(work_pjct.dir_handle, f"{epoch_id}-{branch_name}")

                this_iter = Branch(
                    DirHandle(branch_path, make=True),
                    branch_name=branch_name,
                    chainIDs=old_iter.chainIDs,
                    resnames=branch_resnames,
                    resSeqs=old_iter.resSeqs,
                    parent=old_iter,
                    mutation=mutation)

                init_wt = Path(branch_path, "init_wt.pdb")
                sh.copy(old_pdb, init_wt)
                log.info(f"New mutation: {mutation} on Epoch-Branch: {epoch_id}-{branch_name}")

                # Mutate the PDB
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    try:
                        overlapped_pdb = mutator.on_pdb(
                            PDBStructure.from_path(init_wt),
                            branch_path,
                            mutation=mutation,
                            selection_complex=old_iter.complex.top.selection_complex,
                            selection_wations=old_iter.complex.top.selection_not_complex)
                    except AssertionError as e:
                        # Mutator failed. This position will still be memorized.
                        log.info(f"Mutation of {branch_name} failed. Will try with another one. Backing-up {branch_path} .")
                        print(e, file=sys.stderr)
                        failed_branch_path = Path(work_pjct.dir_handle, f"failed_{epoch_id}-{branch_name}")
                        sh.move(branch_path, failed_branch_path)
                        continue
                remove_overlapping_solvent(
                    overlapped_pdb,
                    mutation.resSeq,
                    Path(branch_path, f"{name}.pdb"),
                    log,
                    use_tleap=work_pjct.config["md"]["use_tleap"],
                )

                # Copy tleap files, if necessary
                work_pjct.get_tleap_into_iter(Path(this_iter.dir_handle))

                this_iter.complex = GROComplex.from_pdb(
                    name=name,
                    input_dir=branch_path,
                    target_chains=work_pjct.config["target"]["chainID"],
                    binder_chains=work_pjct.config["binder"]["chainID"],
                    md_config=work_pjct.config["md"],
                    add_ions=True,
                )
                current_epoch[branch_name] = this_iter
                successful_mutations += 1

            # TODO: check if this works with mutations on different positions
            memorize_mutations(work_pjct, current_epoch, mutations)
    work_pjct.new_epoch(current_epoch)
