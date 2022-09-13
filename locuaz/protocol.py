from pathlib import Path
from projectutils import WorkProject, Epoch, Iteration
from fileutils import DirHandle
from mutationgenerators import mutation_generators
from mutators import mutators
from molecules import split_solute_and_solvent, catenate_pdbs, GROComplex
from gromacsutils import remove_overlapping_waters


def initialize_new_epoch(work_pjct: WorkProject) -> None:
    """initialize_new_epoch(): This is a specific protocol, others will be added

    Args:
        work_pjct (WorkProject): work project
    """
    old_epoch = work_pjct.epochs[-1]
    epoch_id = old_epoch.id + 1
    current_epoch = Epoch(epoch_id, {})

    # Create required mutator
    mutator = mutators[work_pjct.config["protocol"]["mutator"]](
        work_pjct.config["paths"]["mutator"]
    )
    # Create required mutation generator and generate mutation.
    mutation_generator = mutation_generators[work_pjct.config["protocol"]["generator"]](
        old_epoch, work_pjct.config["protocol"]["max_branches"]
    )

    for old_iter_name, mutations in mutation_generator.items():
        old_iter = old_epoch[old_iter_name]
        nonwat_pdb, wation_pdb = split_solute_and_solvent(old_iter.complex)

        for mutation in mutations:
            iter_name, iter_resnames = mutation.new_name_resname(old_iter)
            iter_path = Path(work_pjct.dir_handle, f"{epoch_id}-{iter_name}")
            # Path(str(epoch_id) + "-" + iter_name)
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
            catenate_pdbs(dry_mut_pdb, wation_pdb, pdb_out_path=mut_pdb_fn)
            # Remove the temporary mutated complex that lacks the solvent
            dry_mut_pdb.unlink()

            overlapped_cpx = GROComplex.from_pdb(
                name=over_name,
                input_dir=iter_path,
                target_chains=work_pjct.config["target"]["chainID"],
                binder_chains=work_pjct.config["binder"]["chainID"],
                gmx_bin=work_pjct.config["md"]["gmx_bin"],
            )
            #
            this_iter.complex = remove_overlapping_waters(
                work_pjct.config, overlapped_cpx, mutation.resSeq
            )

            current_epoch[iter_name] = this_iter

    # for top_iter_name in old_epoch.top_iterations:

    #     old_iter = old_epoch[top_iter_name]
    #     mutations = generate_mutations(
    #         old_iter, max_branches=work_pjct.config["protocol"]["max_branches"]
    #     )

    #     nonwat_pdb, wation_pdb = split_solute_and_solvent(old_iter.complex)

    #     for mutation in mutations:
    #         iter_name, iter_resnames = mutation.new_name_resname(old_iter)
    #         iter_path = Path(
    #             str(work_pjct.dir_handle), Path(str(old_epoch.id + 1) + "-" + iter_name)
    #         )
    #         this_iter = Iteration(
    #             DirHandle(iter_path, make=True),
    #             iter_name=iter_name,
    #             chainIDs=old_iter.chainIDs,
    #             resnames=iter_resnames,
    #             resSeqs=old_iter.resSeqs,
    #         )

    #         # Mutate the complex
    #         dry_mut_pdb = mutator(nonwat_pdb, mutation)
    #         # Rejoin the mutated complex with water and ions
    #         over_name = "overlapped_" + work_pjct.config["main"]["name"]
    #         mut_pdb_fn = iter_path / (over_name + ".pdb")
    #         catenate_pdbs(dry_mut_pdb, wation_pdb, pdb_out_path=mut_pdb_fn)
    #         # Remove the temporary mutated complex that lacks the solvent
    #         dry_mut_pdb.unlink()

    #         overlapped_cpx = GROComplex.from_pdb(
    #             name=over_name,
    #             input_dir=iter_path,
    #             target_chains=work_pjct.config["target"]["chainID"],
    #             binder_chains=work_pjct.config["binder"]["chainID"],
    #             gmx_bin=work_pjct.config["md"]["gmx_bin"],
    #         )
    #         #
    #         this_iter.complex = remove_overlapping_waters(
    #             work_pjct.config, overlapped_cpx, mutation.resSeq
    #         )

    #         current_epoch[iter_name] = this_iter
    work_pjct.new_epoch(current_epoch)
