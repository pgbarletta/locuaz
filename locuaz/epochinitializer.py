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
from mutationcreation import MutationCreator
from locuaz.basemutator import memorize_mutations
from locuaz.mutators import mutators
from locuaz.projectutils import WorkProject, Epoch, Branch


def initialize_new_epoch(work_pjct: WorkProject, log: Logger) -> None:
    """initialize_new_epoch(): This is a specific protocol, others will be added

    Args:
        work_pjct (WorkProject):
        log (Logger):
    """
    config = work_pjct.config
    name = config["main"]["name"]
    old_epoch = work_pjct.epochs[-1]
    epoch_id = old_epoch.id + 1
    current_epoch = Epoch(epoch_id, branches={}, nvt_done=False, npt_done=False)

    mutator = mutators[config["mutation"]["mutator"]](
        config["paths"]["mutator"], config["mutation"]["reconstruct_radius"]
    )
    # TODO: Deprecate
    if config.get("generation"):
        generator = mutation_generators[config["generation"]["generator"]]
    successful_mutations = 0

    if config["protocol"]["constant_width"]:
        branches = config["protocol"]["new_branches"]
    else:
        branches = config["protocol"]["new_branches"] * len(old_epoch.top_branches)

    # Usually, this `while` would only be executed once, unless the Mutator program fails to perform a mutation.
    while successful_mutations < branches:
        # TODO: Deprecate
        if config.get("creation"):
            mutation_generator_creator = MutationCreator(old_epoch,
                                                         branches - successful_mutations,
                                                         config["creation"],
                                                         excluded_aas=work_pjct.get_mem_aminoacids(),
                                                         excluded_pos=work_pjct.get_mem_positions(),
                                                         use_tleap=config["md"]["use_tleap"])
        else:
            # TODO: Deprecate
            mutation_generator_creator = generator(
                old_epoch,
                branches - successful_mutations,
                excluded_aas=work_pjct.get_mem_aminoacids(),
                excluded_pos=work_pjct.get_mem_positions(),
                use_tleap=config["md"]["use_tleap"],
                logger=log,
                probe_radius=config["generation"]["probe_radius"])

        for old_branch_name, mutations in mutation_generator_creator.items():
            old_branch = old_epoch.top_branches[old_branch_name]

            if config["md"]["use_tleap"]:
                # GROMACS renumbers resSeqs to strided numbering. If using Amber's continuous
                # numbering, this will result in the wrong mutating_resSeq.
                # Backup the PDB before running pdb4amber
                pdb_path = Path(old_branch.complex.pdb)
                pre_fix_pdb = Path(old_branch.dir_handle, f"preAmberPDBFixer_{pdb_path.stem}.pdb")
                sh.move(pdb_path, pre_fix_pdb)

                old_pdb = fix_pdb(pre_fix_pdb, pdb_path)
            else:
                old_pdb = old_branch.complex.pdb

            for mutation in mutations:
                branch_name, branch_resnames = old_branch.generate_name_resname(mutation)
                branch_path = Path(work_pjct.dir_handle, f"{epoch_id}-{branch_name}")

                this_iter = Branch(
                    DirHandle(branch_path, make=True),
                    branch_name=branch_name,
                    chainIDs=old_branch.chainIDs,
                    resnames=branch_resnames,
                    resSeqs=old_branch.resSeqs,
                    parent=old_branch,
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
                            selection_complex=old_branch.complex.top.selection_complex,
                            selection_wations=old_branch.complex.top.selection_not_complex)
                    except AssertionError as e:
                        # Mutator failed. This position will still be memorized.
                        log.info(
                            f"Mutation of {branch_name} failed. Will try with another one. Backing-up {branch_path} .")
                        print(e, file=sys.stderr)
                        failed_branch_path = Path(work_pjct.dir_handle, f"failed_{epoch_id}-{branch_name}")
                        sh.move(branch_path, failed_branch_path)
                        continue
                remove_overlapping_solvent(
                    overlapped_pdb,
                    mutation.resSeq,
                    Path(branch_path, f"{name}.pdb"),
                    log,
                    use_tleap=config["md"]["use_tleap"],
                )

                # Copy tleap files, if necessary
                work_pjct.get_tleap_into_iter(Path(this_iter.dir_handle))

                this_iter.complex = GROComplex.from_pdb(
                    name=name,
                    input_dir=branch_path,
                    target_chains=config["target"]["chainID"],
                    binder_chains=config["binder"]["chainID"],
                    md_config=config["md"],
                    add_ions=True,
                )
                current_epoch[branch_name] = this_iter
                successful_mutations += 1

            # TODO: check if this works with mutations on different positions
            memorize_mutations(work_pjct, current_epoch, mutations)
    work_pjct.new_epoch(current_epoch)
