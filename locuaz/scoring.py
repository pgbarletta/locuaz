import time
from typing import Dict
import logging
from pathlib import Path

from utils_scoring import extract_pdbs
from projectutils import WorkProject
from fileutils import DirHandle
from primitives import launch_biobb
from biobb_analysis.gromacs.gmx_trjconv_str_ens import GMXTrjConvStrEns
from biobb_analysis.gromacs.gmx_image import GMXImage
from scoringfunctions import *


def initialize_scoring_folder(work_pjct: WorkProject, iter_name: str):
    this_iter = work_pjct.epochs[-1][iter_name]
    current_dir = this_iter.dir_handle.dir_path
    this_iter.score_dir = DirHandle(current_dir / "scoring", make=True, replace=True)
    gmx_bin: str = work_pjct.config["md"]["gmx_bin"]

    # First, remove PBC
    pbc_trj = Path(this_iter.score_dir, "pbc_" + this_iter.complex.name + ".xtc")
    remove_box = GMXImage(
        input_traj_path=str(this_iter.complex.tra.file.path),
        input_top_path=str(this_iter.complex.tpr.file.path),
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

    # Then, extract complex PDBs
    # Zip filename with the extracted PDBs
    ens_of_pdbs = Path(
        this_iter.score_dir, "ensemble_" + this_iter.complex.name + ".zip"
    )
    # this only works when there're no lipids or glyco-stuff
    # IDK why this works. The input TPR has waters, but the input traj doesn't
    get_complex = GMXTrjConvStrEns(
        input_traj_path=str(pbc_trj),
        input_top_path=str(this_iter.complex.tpr.file.path),
        output_str_ens_path=str(ens_of_pdbs),
        properties={"gmx_path": gmx_bin, "selection": "Protein"},
    )
    launch_biobb(get_complex)
    nframes_target = extract_pdbs(ens_of_pdbs, "complex")

    # Extract target PDBs
    get_target = GMXTrjConvStrEns(
        input_traj_path=str(pbc_trj),
        input_top_path=str(this_iter.complex.tpr.file.path),
        input_index_path=str(this_iter.complex.ndx.path),
        output_str_ens_path=str(ens_of_pdbs),
        properties={"gmx_path": gmx_bin, "selection": "target"},
    )
    launch_biobb(get_target)
    nframes_complex = extract_pdbs(ens_of_pdbs, "target")

    # Extract binder PDBs
    get_binder = GMXTrjConvStrEns(
        input_traj_path=str(pbc_trj),
        input_top_path=str(this_iter.complex.tpr.file.path),
        input_index_path=str(this_iter.complex.ndx.path),
        output_str_ens_path=str(ens_of_pdbs),
        properties={"gmx_path": gmx_bin, "selection": "binder"},
    )
    launch_biobb(get_binder)
    nframes_binder = extract_pdbs(ens_of_pdbs, "binder")

    assert (
        nframes_target == nframes_binder == nframes_complex
    ), f"Number of frames in trajectory don't match. This shouldn't happen."

    return nframes_complex


def score(work_pjct: WorkProject, iter_name: str, nframes: int) -> None:
    start = time.time()
    this_iter = work_pjct.epochs[-1][iter_name]

    for sf_name, scorer in work_pjct.scorers.items():
        scores = scorer(nframes=nframes, frames_path=this_iter.score_dir)
        if sf_name == "bluues":
            promedio = this_iter.set_score("bluues", scores[0])
            logging.info(f"{sf_name} average score: {promedio}")
            promedio = this_iter.set_score("bmf", scores[1])
            logging.info(f"bmf average score: {promedio}")
        else:
            promedio = this_iter.set_score(sf_name, scores)
            logging.info(f"{sf_name} average score: {promedio}")

    this_iter.write_down_scores()
    logging.info(
        f"Time elapsed during {iter_name}'s {nframes} frames scoring: {time.time() - start}"
    )
