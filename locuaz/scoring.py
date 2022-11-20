import time
from typing import Dict
import logging
from pathlib import Path
import concurrent.futures as cf

from fileutils import DirHandle
from complex import GROComplex
from projectutils import WorkProject, Iteration
from utils_scoring import extract_pdbs, join_target_binder, rm_aux_scoring_files
from primitives import launch_biobb
from gromacsutils import image_traj
from biobb_analysis.gromacs.gmx_trjconv_str_ens import GMXTrjConvStrEns
from scoringfunctions import *


def initialize_scoring_folder(iteration: Iteration, config: Dict) -> int:
    # No amber support for now.
    assert isinstance(iteration.complex, GROComplex)

    iteration.score_dir = DirHandle(Path(iteration, "scoring"), make=True, replace=True)
    gmx_bin: str = config["md"]["gmx_bin"]

    # First, fix all the imaging issues
    fix_trj_fn = Path(iteration.score_dir, "fix_" + iteration.complex.name + ".xtc")
    fix_trj = image_traj(iteration.complex, fix_trj_fn, gmx_bin)

    # Zip filename with the extracted PDBs
    ens_of_pdbs = Path(
        iteration.score_dir, "ensemble_" + iteration.complex.name + ".zip"
    )
    # Target
    get_target = GMXTrjConvStrEns(
        input_traj_path=str(fix_trj),
        input_top_path=str(iteration.complex.tpr.file.path),
        input_index_path=str(iteration.complex.ndx.path),
        output_str_ens_path=str(ens_of_pdbs),
        properties={"gmx_path": gmx_bin, "selection": "target"},
    )
    launch_biobb(get_target)
    nframes = extract_pdbs(ens_of_pdbs, "target", new_chainID="A")

    # Extract binder PDBs
    get_binder = GMXTrjConvStrEns(
        input_traj_path=str(fix_trj),
        input_top_path=str(iteration.complex.tpr.file.path),
        input_index_path=str(iteration.complex.ndx.path),
        output_str_ens_path=str(ens_of_pdbs),
        properties={"gmx_path": gmx_bin, "selection": "binder"},
    )
    launch_biobb(get_binder)
    nframes_binder = extract_pdbs(ens_of_pdbs, "binder", new_chainID="B")

    # Complex PDBs
    assert nframes == nframes_binder
    join_target_binder(Path(iteration.score_dir), nframes)

    return nframes


def score_frames(work_pjct: WorkProject, iteration: Iteration, nframes: int) -> None:
    start = time.time()
    log = logging.getLogger(f"{work_pjct.name}")

    for sf_name, scorer in work_pjct.scorers.items():
        try:
            scores = scorer(nframes=nframes, frames_path=Path(iteration.score_dir))
        except cf.TimeoutError as e:
            raise e

        assert (
            scores is not None
        ), f"This shouldn't happen. Iteration: {iteration.score_dir}"

        if sf_name == "bluues":
            promedio = iteration.set_score("bluues", scores[0])
            log.info(f"{sf_name} average score: {promedio:.3f}")
            promedio = iteration.set_score("bmf", scores[1])
            log.info(f"bmf average score: {promedio:.3f}")
        else:
            promedio = iteration.set_score(sf_name, scores)
            log.info(f"{sf_name} average score: {promedio:.3f}")

    if not work_pjct.config["main"]["debug"]:
        log.info(
            "Removing PDB, PSF, PQR frames and auxiliary scoring files. Set `--debug` to skip this."
        )
        rm_aux_scoring_files(iteration.score_dir, work_pjct.scorers.keys(), nframes)

    iteration.write_down_scores()
    log.info(
        f"Time elapsed during {iteration.iter_name}'s {nframes} frames scoring: {time.time() - start}"
    )


def score(work_pjct: WorkProject, iteration: Iteration) -> None:
    log = logging.getLogger(f"{work_pjct.name}")
    try:
        iteration.read_scores(work_pjct.scorers.keys(), log)
        log.info("Read old scores.")
    except FileNotFoundError as e:
        log.info("Splitting NPT trajectory in frames.")
        nframes = initialize_scoring_folder(iteration, work_pjct.config)
        score_frames(work_pjct, iteration, nframes)
