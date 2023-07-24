import concurrent.futures as cf
import logging
import sys
import time
from pathlib import Path
from typing import Optional, Tuple

from biobb_analysis.gromacs.gmx_trjconv_str_ens import GMXTrjConvStrEns

from locuaz.fileutils import DirHandle
from locuaz.gromacsutils import image_traj
from locuaz.primitives import launch_biobb
from locuaz.projectutils import WorkProject, Branch
from locuaz.utils_scoring import extract_pdbs, join_target_binder, rm_aux_scoring_files
from locuaz.statistics import run_stats


def initialize_scoring_folder(
        branch: Branch, config: dict, *, log: Optional[logging.Logger] = None
) -> Tuple[int, int]:
    """
    Creates scoring folder inside the branch dir, fixes PBC issues with the original NPT trajectory
    to create a 'fix_{name}.xtc' trajectory and a 'fix_{name}.pdb'.
    'target' and 'binder' will get chainIDs of 'A' and 'B', no matter  their original chainIDs
    or number of chains. This is to prevent scorers from choking.
    Args:
        branch (Branch): Branch object
        config (dict): input config file
        log (logging.Logger): logger

    Returns:
        nframes(int): number of frames
    """
    log.info("Splitting NPT trajectory in frames.")
    branch.score_dir = DirHandle(Path(branch, "scoring"), make=True, replace=True)
    gmx_bin: str = "gmx"

    # First, fix all the imaging issues
    fix_trj_fn = Path(branch.score_dir, "fix_" + branch.complex.name + ".xtc")
    fix_trj = image_traj(branch.complex, fix_trj_fn, config["md"]["use_tleap"], gmx_bin)

    # Zip filename with the extracted PDBs
    ens_of_pdbs = Path(
        branch.score_dir, "ensemble_" + branch.complex.name + ".zip"
    )
    # Target
    get_target = GMXTrjConvStrEns(
        input_traj_path=str(fix_trj),
        input_top_path=str(branch.complex.tpr.file.path),
        input_index_path=str(branch.complex.ndx.path),
        output_str_ens_path=str(ens_of_pdbs),
        properties={"binary_path": gmx_bin, "selection": "target"},
    )
    launch_biobb(get_target)
    start, end = extract_pdbs(
        ens_of_pdbs,
        "target",
        nprocs=config["scoring"]["nthreads"],
        new_chainID="A",
        start=config["scoring"]["start"],
        end=config["scoring"]["end"],
        allowed_nonstandard_residues=set(config["scoring"]["allowed_nonstandard_residues"]),
        log=log,
    )

    # Extract binder PDBs
    get_binder = GMXTrjConvStrEns(
        input_traj_path=str(fix_trj),
        input_top_path=str(branch.complex.tpr.file.path),
        input_index_path=str(branch.complex.ndx.path),
        output_str_ens_path=str(ens_of_pdbs),
        properties={"binary_path": gmx_bin, "selection": "binder"},
    )
    launch_biobb(get_binder)
    start_binder, end_binder = extract_pdbs(
        ens_of_pdbs,
        "binder",
        nprocs=config["scoring"]["nthreads"],
        new_chainID="B",
        start=config["scoring"]["start"],
        end=config["scoring"]["end"],
        allowed_nonstandard_residues=set(config["scoring"]["allowed_nonstandard_residues"]),
        log=log,
    )

    # Complex PDBs
    assert start == start_binder, f"target start frame ({start}) and binder start frame ({start_binder}) " \
                                  f"don't match. This shouldn't happen."
    assert end == end_binder, f"target end frame ({end}) and binder end frame ({end_binder}) " \
                              f"don't match. This shouldn't happen."
    join_target_binder(Path(branch.score_dir), start=start, end=end)

    return start, end


def score_frames(work_pjct: WorkProject, branch: Branch, *, start: int, end: int) -> None:
    log = logging.getLogger(f"{work_pjct.name}")
    log.info(f"Beginning scoring.")
    for sf_name, scorer in work_pjct.scorers.items():
        try:
            scores = scorer(start=start, end=end, frames_path=Path(branch.score_dir), cpx=branch.complex)
        except cf.TimeoutError as e:
            raise e

        assert scores is not None, f"This shouldn't happen. Branch: {branch.score_dir}"

        avg_val = branch.set_score(sf_name, scores)
        log.info(f"{sf_name} average score: {avg_val:.3f}")

    if not work_pjct.config["main"]["debug"]:
        log.info("Removing PDB, PSF, PQR frames and auxiliary scoring files. Set `--debug` to skip this.")
        rm_aux_scoring_files(branch.score_dir, work_pjct.scorers.keys(), start=start, end=end)

    branch.write_down_scores()


def discard_branch(work_pjct: WorkProject, branch: Branch) -> None:
    log = logging.getLogger(f"{work_pjct.name}")

    for sf_name, _ in work_pjct.scorers.items():
        branch.set_score(sf_name, [sys.maxsize, sys.maxsize])
        log.info(f"{sf_name} nullifying score.")
    # Initialize the scoring folder
    branch.score_dir = DirHandle(Path(branch, "scoring"), make=True, replace=True)
    branch.write_down_scores()


def score(work_pjct: WorkProject, branch: Branch) -> None:
    log = logging.getLogger(f"{work_pjct.name}")

    if branch.read_scores(work_pjct.scorers.keys(), log):
        log.info("Read old scores.")
        return
    else:
        # TODO branch.outside_box should be saved on tracking
        if branch.outside_box:
            # Discard this branch.
            discard_branch(work_pjct, branch)
        else:
            start_time = time.time()
            start, end = initialize_scoring_folder(branch, work_pjct.config, log=log)
            run_stats(branch, work_pjct.config, start=start, end=end, log=log)
            score_frames(work_pjct, branch, start=start, end=end)
            log.info(
                f"Time elapsed during {branch.branch_name}'s {end - start} "
                f"frames scoring: {time.time() - start_time}")
