from pathlib import Path
from projectutils import WorkProject
from biobb_md.gromacs.pdb2gmx import pdb2gmx
from biobb_md.gromacs.solvate import solvate
from molecules import ZipTopology
from biobb_md.gromacs.grompp import grompp
from biobb_md.gromacs.grompp_mdrun import grompp_mdrun
from biobb_md.gromacs.mdrun import mdrun
from molecules import (
    GROComplex,
    ZipTopology,
    GROStructure,
    PDBStructure,
    XtcTrajectory,
    TPRFile,
)
from typing import Dict
from biobb_analysis.gromacs.gmx_trjconv_str import gmx_trjconv_str
from runutils import MDrun


def prepare_for_mdrun(work_dir: WorkProject, iter_name: str, maxsol: int = 1):
    if maxsol == 0:
        print(
            "prepare_for_mdrun(): Warning, maxsol == 0, gmx solvate "
            "we'll add as many waters as it wants."
        )

    config = work_dir.config
    this_iter = work_dir.epochs[-1][iter_name]
    system_name = this_iter.complex.pdb.file.name
    current_dir = this_iter.dir_handle.dir_path

    # PDB2GMX
    dry_gro = str(current_dir / ("dry_" + system_name + ".gro"))
    dry_top_zip = str(current_dir / "dry_topol.zip")
    props = {
        "gmx_path": str(config["md"]["gmx_bin"]),
        "water_type": "tip3p",
        "force_field": "amber99sb-ildn",
        "ignh": True,
    }
    pdb2gmx(
        input_pdb_path=str(this_iter.complex.pdb.file.path),
        output_gro_path=dry_gro,
        output_top_zip_path=dry_top_zip,
        properties=props,
    )

    # SOLVATE
    wet_gro = str(current_dir / ("wet_" + system_name + ".gro"))
    wet_top_zip = str(current_dir / "wet_topol.zip")
    props = {"gmx_path": str(config["md"]["gmx_bin"]), "dev": f"-maxsol {maxsol}"}
    solvate(
        input_solute_gro_path=dry_gro,
        output_gro_path=wet_gro,
        input_top_zip_path=dry_top_zip,
        output_top_zip_path=wet_top_zip,
        properties=props,
    )

    a = ZipTopology.from_path_with_chains(
        wet_top_zip,
        target_chains=work_dir.config["target"]["chainID"],
        binder_chains=work_dir.config["binder"]["chainID"],
    )
    b = GROStructure.from_path(wet_gro)
    c = PDBStructure.from_path(this_iter.complex.pdb.file.path)
    this_iter.complex = GROComplex(system_name, current_dir, c, a, b)
    this_iter.complex.update_all_ndxs()


def run_min_nvt_npt(
    work_dir: WorkProject, iter_name: str, config: Dict, gpu_id: int = 0
):

    this_iter = work_dir.epochs[-1][iter_name]
    system_name = this_iter.complex.pdb.file.name
    current_dir = this_iter.dir_handle.dir_path

    # MDRUN: MIN
    min_tpr = str(current_dir / ("min_" + system_name + ".tpr"))
    grompp(
        input_mdp_path=str(this_iter.mdps["min_mdp"].path),
        input_gro_path=str(this_iter.complex.gro.file.path),
        input_top_zip_path=str(this_iter.complex.top.file.path),
        output_tpr_path=str(min_tpr),
        properties={"gmx_path": str(config["md"]["gmx_bin"])},
    )

    min_trr = str(current_dir / ("min_" + system_name + ".trr"))
    min_gro = str(current_dir / ("min_" + system_name + ".gro"))
    min_edr = str(current_dir / ("min_" + system_name + ".edr"))
    min_log = str(current_dir / ("min_" + system_name + ".log"))

    props = {
        "gmx_path": str(config["md"]["gmx_bin"]),
        "num_threads_omp": config["md"]["omp_procs"],
        "num_threads_mpi": 1,
        # "use_gpu": True,
        "dev": "-nb gpu -pin on -pinoffset 0 -pinstride 1",
    }

    mdrun(
        input_tpr_path=min_tpr,
        output_trr_path=min_trr,
        output_gro_path=min_gro,
        output_edr_path=min_edr,
        output_log_path=min_log,
        properties=props,
    )

    # MDRUN: NVT
    nvt_trr = str(current_dir / ("nvt_" + system_name + ".trr"))
    nvt_xtc = str(current_dir / ("nvt_" + system_name + ".xtc"))
    nvt_gro = str(current_dir / ("nvt_" + system_name + ".gro"))
    nvt_edr = str(current_dir / ("nvt_" + system_name + ".edr"))
    nvt_log = str(current_dir / ("nvt_" + system_name + ".log"))
    props = {
        "gmx_path": str(config["md"]["gmx_bin"]),
        "num_threads_omp": config["md"]["omp_procs"],
        "num_threads_mpi": config["md"]["mpi_procs"],
        "use_gpu": True,
        "gpu_id": gpu_id,
        "dev": "-bonded gpu -pin on -pinoffset 0 -pinstride 1 -pmefft gpu",
    }

    grompp_mdrun(
        input_gro_path=min_gro,
        input_top_zip_path=str(this_iter.complex.top.file.path),
        input_mdp_path=str(this_iter.mdps["nvt_mdp"].path),
        output_trr_path=nvt_trr,
        output_xtc_path=nvt_xtc,
        output_gro_path=nvt_gro,
        output_edr_path=nvt_edr,
        output_log_path=nvt_log,
        properties=props,
    )

    # MDRUN: NPT
    npt_tpr = str(current_dir / ("npt_" + system_name + ".tpr"))
    grompp(
        input_mdp_path=str(this_iter.mdps["npt_mdp"].path),
        input_gro_path=nvt_gro,
        input_top_zip_path=str(this_iter.complex.top.file.path),
        output_tpr_path=str(npt_tpr),
        properties={},
    )

    npt_trr = str(current_dir / ("npt_" + system_name + ".trr"))
    npt_xtc = str(current_dir / ("npt_" + system_name + ".xtc"))
    npt_gro = str(current_dir / ("npt_" + system_name + ".gro"))
    npt_edr = str(current_dir / ("npt_" + system_name + ".edr"))
    npt_log = str(current_dir / ("npt_" + system_name + ".log"))
    props = {
        "gmx_path": str(config["md"]["gmx_bin"]),
        "num_threads_omp": config["md"]["omp_procs"],
        "num_threads_mpi": config["md"]["mpi_procs"],
        "use_gpu": True,
        "gpu_id": gpu_id,
        "dev": "-bonded gpu -pin on -pinoffset 0 -pinstride 1 -pmefft gpu",
    }

    mdrun(
        input_tpr_path=npt_tpr,
        output_trr_path=npt_trr,
        output_xtc_path=npt_xtc,
        output_gro_path=npt_gro,
        output_edr_path=npt_edr,
        output_log_path=npt_log,
        properties=props,
    )

    # Get last frame in PDB:
    npt_pdb = str(current_dir / ("npt_" + system_name + ".pdb"))
    props = {"gmx_path": str(config["md"]["gmx_bin"]), "selection": "System"}
    gmx_trjconv_str(
        input_structure_path=npt_gro,
        input_top_path=npt_tpr,
        output_str_path=npt_pdb,
        properties=props,
    )

    this_iter.complex = GROComplex(
        system_name,
        current_dir,
        pdb=PDBStructure.from_path(npt_pdb),
        top=this_iter.complex.top,
        gro=GROStructure.from_path(Path(npt_gro)),
    )

    this_iter.complex.tra = XtcTrajectory.from_path(npt_xtc)
    this_iter.complex.tpr = TPRFile.from_path(npt_tpr)
    this_iter.complex.update_all_ndxs()
