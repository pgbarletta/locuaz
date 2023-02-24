import shutil as sh
import warnings
from logging import Logger
from pathlib import Path

from amberutils import fix_pdb
from complex import GROComplex
from fileutils import DirHandle
from gromacsutils import remove_overlapping_solvent
from mutationgenerators import mutation_generators
from mutator import memorize_mutations
from mutators import mutators
from projectutils import WorkProject, Epoch, Iteration


def initialize_new_epoch(work_pjct: WorkProject, log: Logger) -> None:
    """initialize_new_epoch(): This is a specific protocol, others will be added

    Args:
        work_pjct (WorkProject):
        log (Logger):
    """
    name = work_pjct.config["main"]["name"]
    old_epoch = work_pjct.epochs[-1]
    epoch_id = old_epoch.id + 1
    current_epoch = Epoch(epoch_id, iterations={}, nvt_done=False, npt_done=False)

    # Create required mutator
    mutator = mutators[work_pjct.config["protocol"]["mutator"]](
        work_pjct.config["paths"]["mutator"]
    )
    # Create required mutation generator and generate mutation.
    generator = mutation_generators[work_pjct.config["protocol"]["generator"]]
    mutation_generator = generator(
        old_epoch,
        work_pjct.config["protocol"]["branches"],
        excluded_aas=work_pjct.get_mem_aminoacids(),
        excluded_pos=work_pjct.get_mem_positions(),
        use_tleap=work_pjct.config["md"]["use_tleap"],
        logger=log
    )

    for old_iter_name, mutations in mutation_generator.items():
        old_iter = old_epoch.top_iterations[old_iter_name]

        # GROMACS renumbers resSeqs to strided numbering. If using Amber's continuous
        # numbering, this will result in the wrong mutating_resSeq.
        if work_pjct.config["md"]["use_tleap"]:
            # Backup the PDB before runing pdb4amber
            pdb_path = Path(old_iter.complex.pdb)
            pre_fix_pdb = Path(
                old_iter.dir_handle, f"preAmberPDBFixer_{pdb_path.stem}.pdb"
            )
            sh.move(pdb_path, pre_fix_pdb)

            old_pdb = fix_pdb(pre_fix_pdb, pdb_path)
        else:
            old_pdb = old_iter.complex.pdb

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
            log.info(
                f"New mutation: {mutation} on Epoch-Iteration: {epoch_id}-{iter_name}"
            )

            # Mutate the PDB
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                overlapped_pdb = mutator.on_pdb(
                    old_pdb,
                    iter_path,
                    mutation=mutation,
                    selection_protein=old_iter.complex.top.selection_protein,
                    selection_wations=old_iter.complex.top.selection_not_protein,
                )
            remove_overlapping_solvent(
                overlapped_pdb,
                mutation.resSeq,
                Path(iter_path, f"{name}.pdb"),
                log,
                use_tleap=work_pjct.config["md"]["use_tleap"],
            )

            # Copy tleap files, if necessary
            work_pjct.get_tleap_into_iter(Path(this_iter.dir_handle))

            this_iter.complex = GROComplex.from_pdb(
                name=name,
                input_dir=iter_path,
                target_chains=work_pjct.config["target"]["chainID"],
                binder_chains=work_pjct.config["binder"]["chainID"],
                md_config=work_pjct.config["md"],
                add_ions=True,
            )
            current_epoch[iter_name] = this_iter

        # TODO: check if this works with mutations on different positions
        memorize_mutations(work_pjct, current_epoch, mutations)
    work_pjct.new_epoch(current_epoch)
