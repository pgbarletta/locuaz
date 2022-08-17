import time
from typing import Dict
from utils_scoring import extract_pdbs
from projectutils import WorkProject
from fileutils import DirHandle
from biobb_analysis.gromacs.gmx_trjconv_str_ens import gmx_trjconv_str_ens
from scoringfunctions import *


def initialize_scoring_folder(work_pjct: WorkProject, iter_name: str, config: Dict):
    print("--- Creating scoring folder and extracting PDBs from NPT trajectory ---")

    this_iter = work_pjct.epochs[-1][iter_name]
    current_dir = this_iter.dir_handle.dir_path
    this_iter.score_dir = DirHandle(current_dir / "scoring", make=True, replace=True)

    # Extract target PDBs
    ens_of_pdbs = this_iter.score_dir.dir_path / (this_iter.complex.name + ".zip")
    props = {"selection": "target"}
    gmx_trjconv_str_ens(
        input_traj_path=str(this_iter.complex.tra.file.path),
        input_top_path=str(this_iter.complex.tpr.file.path),
        input_index_path=str(this_iter.complex.complex_ndx.path),
        output_str_ens_path=str(ens_of_pdbs),
        properties=props,
    )
    nframes_target = extract_pdbs(ens_of_pdbs, "target")

    # Extract binder PDBs
    ens_of_pdbs = this_iter.score_dir.dir_path / (this_iter.complex.name + ".zip")
    props = {"selection": "binder"}
    gmx_trjconv_str_ens(
        input_traj_path=str(this_iter.complex.tra.file.path),
        input_top_path=str(this_iter.complex.tpr.file.path),
        input_index_path=str(this_iter.complex.complex_ndx.path),
        output_str_ens_path=str(ens_of_pdbs),
        properties=props,
    )
    nframes_binder = extract_pdbs(ens_of_pdbs, "binder")

    # Extract complex PDBs
    ens_of_pdbs = this_iter.score_dir.dir_path / (this_iter.complex.name + ".zip")
    props = {"selection": "Protein"}
    gmx_trjconv_str_ens(
        input_traj_path=str(this_iter.complex.tra.file.path),
        input_top_path=str(this_iter.complex.tpr.file.path),
        # input_index_path=str(this_iter.complex.complex_ndx.path),
        output_str_ens_path=str(ens_of_pdbs),
        properties=props,
    )
    nframes_complex = extract_pdbs(ens_of_pdbs, "complex")

    assert (
        nframes_target == nframes_binder == nframes_complex
    ), f"Number of frames in trajectory don't match. This shouldn't happen."

    return nframes_complex


def score(work_pjct: WorkProject, iter_name: str, config: Dict, nframes: int) -> None:
    start = time.time()
    this_iter = work_pjct.epochs[-1][iter_name]

    for sf_name, scorer in work_pjct.scorers.items():
        scores = scorer(nframes=nframes, frames_path=this_iter.score_dir)
        if sf_name == "bluues_scorer":
            promedio = this_iter.set_score("bluues", scores[0])
            print(f"bluues: {promedio}")
            promedio = this_iter.set_score("bmf", scores[1])
            print(f"bmf: {promedio}")
        else:
            promedio = this_iter.set_score(sf_name, scores)
            print(promedio)
        print(f"-------")

    this_iter.write_down_scores()
    end = time.time()
    print(end - start)
