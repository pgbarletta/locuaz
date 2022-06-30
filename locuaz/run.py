from pathlib import Path
from attrs import define, field
from projectutils import WorkProject  # type: ignore
from biobb_md.gromacs.pdb2gmx import pdb2gmx  # type: ignore
from biobb_md.gromacs.solvate import solvate  # type: ignore
from molecules import ZipTopology  # type: ignore
from biobb_md.gromacs.grompp import grompp  # type: ignore
from biobb_md.gromacs.genion import genion  # type: ignore
from biobb_md.gromacs.grompp_mdrun import grompp_mdrun  # type: ignore


# from runutils import RunOptions


@define
class RunOptions:
    nprocs: int = field(converter=int, default=4, kw_only=True)
    ngpus: int = field(converter=int, default=1, kw_only=True)


def prepare_for_mdrun(work_dir: WorkProject, iter_name: str, maxsol: int = 1):
    if maxsol == 0:
        print(
            "run_min_nvt_npt(): Warning, maxsol == 0, gmx solvate "
            "we'll add as many waters as it wants."
        )

    this_iter = work_dir.epochs[-1][iter_name]
    system_name = this_iter.complex.pdb.file.name
    current_dir = this_iter.dir_handle.dir_path

    # PDB2GMX
    dry_gro = str(current_dir / ("dry_" + system_name + ".gro"))
    dry_top_zip = str(current_dir / "dry_topol.zip")
    props = {"water_type": "tip3p", "force_field": "amber99sb-ildn", "ignh": True}
    pdb2gmx(
        input_pdb_path=str(this_iter.complex.pdb.file.path),
        output_gro_path=dry_gro,
        output_top_zip_path=dry_top_zip,
        properties=props,
    )

    # SOLVATE
    wet_gro = str(current_dir / ("wet_" + system_name + ".gro"))
    wet_top_zip = str(current_dir / "wet_topol.zip")
    props = {"dev": f"-maxsol {maxsol}"}
    solvate(
        input_solute_gro_path=dry_gro,
        output_gro_path=wet_gro,
        input_top_zip_path=dry_top_zip,
        output_top_zip_path=wet_top_zip,
        properties=props,
    )
    # IONS
    # ion_tpr = str(current_dir / ("ion_" + system_name + ".tpr"))
    # props = {"simulation_type": "minimization", "maxwarn": 1}

    # grompp(
    #     input_gro_path=wet_gro,
    #     input_top_zip_path=wet_top_zip,
    #     output_tpr_path=ion_tpr,
    #     properties=props,
    # )

    # genion_gro = str(current_dir / ("gen_" + system_name + ".gro"))
    # genion_top_zip = str(current_dir / "gen_topol.zip")
    # props = {"neutral": True}
    # genion(
    #     input_tpr_path=ion_tpr,
    #     input_top_zip_path=wet_top_zip,
    #     output_gro_path=genion_gro,
    #     output_top_zip_path=genion_top_zip,
    #     properties=props,
    # )


def run_min_nvt_npt(work_dir: WorkProject, iter_name: str, run_options: RunOptions):

    this_iter = work_dir.epochs[-1][iter_name]
    system_name = this_iter.complex.pdb.file.name
    current_dir = this_iter.dir_handle.dir_path

    # MDRUN: minimize
    min_trr = str(current_dir / ("min_" + system_name + ".trr"))
    min_gro = str(current_dir / ("min_" + system_name + ".gro"))
    min_edr = str(current_dir / ("min_" + system_name + ".edr"))
    min_log = str(current_dir / ("min_" + system_name + ".log"))
    props = {
        "num_threads_omp": 4,
        # "num_threads_mpi": 1,
        "gmx_path": "gmx",
        "use_gpu": True,
        "dev": "-pin on -pinoffset 0",
    }

    grompp_mdrun(
        input_gro_path=wet_gro,
        input_top_zip_path=wet_top_zip,
        input_mdp_path=str(this_iter.min_mdp.path),
        output_trr_path=min_trr,
        output_gro_path=min_gro,
        output_edr_path=min_edr,
        output_log_path=min_log,
        properties=props,
    )
