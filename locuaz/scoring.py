import time
from typing import Dict
import logging
from pathlib import Path

from utils_scoring import extract_pdbs, join_target_binder, rm_frames
from projectutils import WorkProject, Iteration
from fileutils import DirHandle
from primitives import launch_biobb
from biobb_analysis.gromacs.gmx_trjconv_str_ens import GMXTrjConvStrEns
from biobb_analysis.gromacs.gmx_image import GMXImage
from scoringfunctions import *


def initialize_scoring_folder(iteration: Iteration, config: Dict) -> int:
    iteration.score_dir = DirHandle(Path(iteration, "scoring"), make=True, replace=True)
    gmx_bin: str = config["md"]["gmx_bin"]

    # First, remove PBC
    pbc_trj = Path(iteration.score_dir, "pbc_" + iteration.complex.name + ".xtc")
    remove_box = GMXImage(
        input_traj_path=str(iteration.complex.tra.file.path),
        input_top_path=str(iteration.complex.tpr.file.path),
        output_traj_path=str(pbc_trj),
        properties={
            "gmx_path": gmx_bin,
            "fit_selection": "Protein",
            "center_selection": "Protein",
            "output_selection": "Protein",
            "pbc": "mol",
        },
    )

    launch_biobb(remove_box)

    # Zip filename with the extracted PDBs
    ens_of_pdbs = Path(
        iteration.score_dir, "ensemble_" + iteration.complex.name + ".zip"
    )
    # Target
    get_target = GMXTrjConvStrEns(
        input_traj_path=str(pbc_trj),
        input_top_path=str(iteration.complex.tpr.file.path),
        input_index_path=str(iteration.complex.ndx.path),
        output_str_ens_path=str(ens_of_pdbs),
        properties={"gmx_path": gmx_bin, "selection": "target"},
    )
    launch_biobb(get_target)
    nframes = extract_pdbs(ens_of_pdbs, "target", new_chainID="A")

    # Extract binder PDBs
    get_binder = GMXTrjConvStrEns(
        input_traj_path=str(pbc_trj),
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

    for sf_name, scorer in work_pjct.scorers.items():
        scores = scorer(nframes=nframes, frames_path=iteration.score_dir)
        if sf_name == "bluues":
            promedio = iteration.set_score("bluues", scores[0])
            logging.info(f"{sf_name} average score: {promedio}")
            promedio = iteration.set_score("bmf", scores[1])
            logging.info(f"bmf average score: {promedio}")
        else:
            promedio = iteration.set_score(sf_name, scores)
            logging.info(f"{sf_name} average score: {promedio}")

    if not work_pjct.config["main"]["debug"]:
        logging.info("Removing PDB frames. Set `--debug` flag to skip this.")
        rm_frames(iteration.score_dir, nframes)

    iteration.write_down_scores()
    logging.info(
        f"Time elapsed during {iteration.iter_name}'s {nframes} frames scoring: {time.time() - start}"
    )


def score(work_pjct: WorkProject, iteration: Iteration) -> None:
    try:
        iteration.read_scores(work_pjct.scorers.keys())
        logging.info("Read old scores.")
    except FileNotFoundError as e:
        nframes = initialize_scoring_folder(iteration, work_pjct.config)
        score_frames(work_pjct, iteration, nframes)
