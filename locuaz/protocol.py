from pathlib import Path
from projectutils import WorkProject, Epoch, Iteration
from fileutils import DirHandle
from mutationgenerators import mutation_generators
from mutators import mutators
from mutator import Mutation, memorize_mutations
from molecules import split_solute_and_solvent, catenate_pdbs, GROComplex
from gromacsutils import remove_overlapping_waters


def initialize_new_epoch(work_pjct: WorkProject) -> None:
    """initialize_new_epoch(): This is a specific protocol, others will be added

    Args:
        work_pjct (WorkProject): work project
    """
    old_epoch = work_pjct.epochs[-1]
    epoch_id = old_epoch.id + 1
    current_epoch = Epoch(epoch_id, iterations={}, nvt_done=False, npt_done=False)

    # Create required mutator
    mutator = mutators[work_pjct.config["protocol"]["mutator"]](
        work_pjct.config["paths"]["mutator"]
    )
    # Create required mutation generator and generate mutation.
    mutation_generator = mutation_generators[work_pjct.config["protocol"]["generator"]](
        old_epoch,
        work_pjct.config["protocol"]["max_branches"],
        excluded_aas=work_pjct.get_mem_aminoacids(),
        excluded_pos=work_pjct.get_mem_positions(),
    )

    for old_iter_name, mutations in mutation_generator.items():
        old_iter = old_epoch.top_iterations[old_iter_name]
        # Get the system's box size after the NPT run, to add it later onto the
        # mutated PDB system. The PDB format has less precision for the box parameters
        # than the GRO format, so there may be a difference in the last digit for the
        # lengths (eg: 12.27215 to 12.27210) and the angles (6.13607 to 6.13605).
        # That's why GROComplex.from_pdb() also uses editconf.
        cryst1_record = old_iter.complex.get_cryst1_record()
        nonwat_pdb, wation_pdb = split_solute_and_solvent(old_iter.complex)

        for mutation in mutations:
            iter_name, iter_resnames = mutation.new_name_resname(old_iter)
            iter_path = Path(work_pjct.dir_handle, f"{epoch_id}-{iter_name}")

            this_iter = Iteration(
                DirHandle(iter_path, make=True),
                iter_name=iter_name,
                chainIDs=old_iter.chainIDs,
                resnames=iter_resnames,
                resSeqs=old_iter.resSeqs,
            )

            # Mutate the complex
            dry_mut_pdb = mutator(nonwat_pdb, mutation)
            # Rejoin the mutated complex with water and ions
            over_name = "overlapped_" + work_pjct.config["main"]["name"]
            mut_pdb_fn = iter_path / (over_name + ".pdb")
            mut_pdb = catenate_pdbs(dry_mut_pdb, wation_pdb, pdb_out_path=mut_pdb_fn)
            mut_pdb.set_cryst1_record(cryst1_record)
            # Remove the temporary mutated complex that lacks the solvent
            dry_mut_pdb.unlink()

            overlapped_cpx = GROComplex.from_pdb(
                name=over_name,
                input_dir=iter_path,
                target_chains=work_pjct.config["target"]["chainID"],
                binder_chains=work_pjct.config["binder"]["chainID"],
                md_config=work_pjct.config["md"],
            )
            #
            this_iter.complex = remove_overlapping_waters(
                work_pjct.config, overlapped_cpx, mutation.resSeq
            )

            current_epoch[iter_name] = this_iter

    memorize_mutations(work_pjct, mutations)
    work_pjct.new_epoch(current_epoch)
