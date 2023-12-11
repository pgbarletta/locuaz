import shutil as sh
import sys
import warnings
from logging import Logger
from os import listdir
from pathlib import Path
from typing import List, Dict, Any, Optional
import concurrent.futures as cf
from multiprocessing import Lock

from locuaz.amberutils import fix_pdb
from locuaz.complex import GROComplex
from locuaz.fileutils import DirHandle
from locuaz.gromacsutils import remove_overlapping_solvent
from locuaz.molecules import PDBStructure
from locuaz.mutationgenerators import mutation_generators
from locuaz.mutationcreation import MutationCreator
from locuaz.basemutator import memorize_mutations, BaseMutator
from locuaz.mutators import mutators
from locuaz.projectutils import WorkProject, Epoch, Branch
from locuaz.mutation import Mutation
from locuaz.primitives import MutationError

lock = Lock()


def initialize_new_epoch(work_pjct: WorkProject, log: Logger) -> Epoch:
    """

    Parameters
    ----------
    work_pjct : WorkProject
        Work Project
    log : Logger
        Logger
    Returns
    -------
    None
    """
    config: Dict[str, Any] = work_pjct.config
    old_epoch = work_pjct.epochs[-1]
    new_epoch = Epoch(old_epoch.id + 1, branches={})

    mutator = mutators[config["mutation"]["mutator"]](
        config["paths"]["mutator"],
        config["mutation"]["reconstruct_radius"],
        config["mutation"]["allowed_nonstandard_residues"],
    )

    if config.get("creation"):
        new_branches = config["protocol"]["new_branches"] * config["creation"]["sites"]
    else:
        # TODO: Deprecate
        generator = mutation_generators[config["generation"]["generator"]]
        new_branches = config["protocol"]["new_branches"]

    if not config["protocol"]["constant_width"]:
        new_branches *= len(old_epoch.top_branches)

    # Control variables for the `while` loop.
    successful_mutations = 0
    failed_already = False
    # Will use this to memorize mutations performed on the epoch
    all_mutations: List[Mutation] = []
    # Usually, this `while` would only be executed once, unless the Mutator
    # program fails to perform a mutation or two or more mutations end up generating
    # the same branch.
    while True:
        new_branches -= successful_mutations
        if config.get("creation"):
            mutation_generator_creator = MutationCreator(
                old_epoch.top_branches,
                new_branches,
                config["creation"],
                excluded_sites=work_pjct.get_mem_positions(),
                amber_numbering=config["md"]["use_tleap"],
                logger=log,
            )
        else:
            # TODO: Deprecate
            mutation_generator_creator = generator(
                old_epoch,
                new_branches,
                excluded_aas=work_pjct.get_mem_aminoacids(),
                excluded_pos=work_pjct.get_mem_positions(),
                use_tleap=config["md"]["use_tleap"],
                logger=log,
                probe_radius=config["generation"]["probe_radius"],
            )

        actual_new_branches = sum(
            [len(muts) for muts in mutation_generator_creator.values()]
        )
        with cf.ProcessPoolExecutor(max_workers=config["scoring"]["nthreads"]) as ex:
            futuros_branches = {}
            for old_branch_name, mutations in mutation_generator_creator.items():
                all_mutations.extend(mutations)
                old_branch = old_epoch.top_branches[old_branch_name]
                for mutation in mutations:
                    futu_branch = ex.submit(
                        create_branch,
                        work_pjct.name,
                        old_branch,
                        mutator=mutator,
                        mutation=mutation,
                        md_config=config["md"],
                        tleap_dir=work_pjct.tleap_dir,
                        log=log,
                    )
                    futuros_branches[futu_branch] = (old_branch.branch_name, mutation)
            for futu_branch in cf.as_completed(futuros_branches):
                old_branch_name, mutation = futuros_branches[futu_branch]
                try:
                    branch = futu_branch.result()
                    new_epoch[branch.branch_name] = branch
                    successful_mutations += 1
                except MutationError:
                    log.warning(
                        f"Mutator failed to mutate {old_branch_name} with {mutation=}"
                        "\nWill try again with another mutation on another branch."
                    )
                    continue
                except FileExistsError:
                    log.warning(
                        f"Mutator tried to mutate {old_branch_name} with {mutation=} but "
                        "this would generate a branch that's identical to a recently "
                        " created one."
                    )
                    continue
                except Exception as e:
                    raise e
            # Using `<=` instead of `==` because if the `while` loops runs twice,
            # `actual_new_branches` will probably be much lower than `successful_mutations`
            # since this last variable accumulates the successful mutations from the
            # 2 loop runs.
            if actual_new_branches <= successful_mutations:
                log.info(
                    f"{actual_new_branches} out of {new_branches} branches created."
                )
                break
            else:
                log.info(
                    f"Tried to generate {new_branches} new branches. Generated "
                    f"{actual_new_branches} new mutations, but only generated "
                    f"{successful_mutations} branches."
                )
                if failed_already:
                    log.info("Will move ahead.")
                    break
                else:
                    log.info(
                        f"Will try to generate the {new_branches - successful_mutations} "
                        "missing branches."
                    )
                    failed_already = True
    # TODO: check if this works with mutations on different positions
    memorize_mutations(work_pjct, new_epoch, all_mutations)
    return new_epoch


def create_branch(
    name: str,
    old_branch: Branch,
    *,
    mutator: BaseMutator,
    mutation: Mutation,
    md_config: Dict[str, Any],
    tleap_dir: Optional[DirHandle],
    log: Logger,
) -> Branch:
    branch_name, branch_resnames = old_branch.generate_name_resname(mutation)
    epoch_id = old_branch.epoch_id + 1
    work_dir = Path(old_branch.dir_handle).parent
    branch_path = Path(work_dir, f"{epoch_id}-{branch_name}")

    if md_config["use_tleap"]:
        with lock:
            old_pdb = fix_old_pdb_numbering(old_branch)
    else:
        old_pdb = old_branch.complex.pdb

    new_branch = Branch(
        DirHandle(branch_path, make=True),
        branch_name=branch_name,
        chainIDs=old_branch.chainIDs,
        resnames=branch_resnames,
        resSeqs=old_branch.resSeqs,
        parent=old_branch,
        mutation=mutation,
    )

    init_wt = Path(branch_path, "init_wt.pdb")
    with lock:
        sh.copy(old_pdb, init_wt)
    log.info(
        f"{mutation} on: {old_branch.epoch_id}-{old_branch.branch_name} "
        f"to generate: {epoch_id}-{branch_name}"
    )

    # Mutate the PDB
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            overlapped_pdb = mutator.on_pdb(
                PDBStructure.from_path(init_wt),
                branch_path,
                mutation=mutation,
                selection_complex=old_branch.complex.top.selection_complex,
                selection_wations=old_branch.complex.top.selection_not_complex,
            )
        except AssertionError as e:
            # Mutator failed. This position will still be memorized.
            log.info(
                f"Mutation of {branch_name} failed. Will try with another one. Backing-up {branch_path} ."
            )
            print(e, file=sys.stderr)
            failed_branch_path = Path(work_dir, f"failed_{epoch_id}-{branch_name}")
            sh.move(branch_path, failed_branch_path)
            raise MutationError from e

    remove_overlapping_solvent(
        overlapped_pdb,
        mutation.resSeq,
        Path(branch_path, f"{name}.pdb"),
        log,
        use_tleap=md_config["use_tleap"],
    )

    # Copy tleap files, if necessary
    if md_config["use_tleap"]:
        for file in listdir(Path(tleap_dir)):
            sh.copy(Path(tleap_dir, file), Path(new_branch.dir_handle))

    new_branch.complex = GROComplex.from_pdb(
        name=name,
        input_dir=branch_path,
        target_chains=old_branch.complex.top.target_chains,
        binder_chains=old_branch.complex.top.binder_chains,
        md_config=md_config,
        add_ions=True,
    )
    return new_branch


def fix_old_pdb_numbering(old_branch: Branch) -> PDBStructure:
    """
    GROMACS renumbers resSeqs to strided numbering. If using Amber's
    continuous numbering, this will result in the wrong mutating_resSeq.

    Parameters
    ----------
    old_branch : Branch
        branch about to be mutated

    Returns
    -------
    fixed_old_pdb : PDBStructure
        The Branch's PDB with amber's resSeq numbering scheme.
    """
    # Backup the PDB before running pdb4amber
    pdb_path = Path(old_branch.complex.pdb)
    pre_fix_pdb = Path(old_branch.dir_handle, f"preAmberPDBFixer_{pdb_path.stem}.pdb")
    sh.move(pdb_path, pre_fix_pdb)

    return fix_pdb(pre_fix_pdb, pdb_path)
