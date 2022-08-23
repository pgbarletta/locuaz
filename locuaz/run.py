from projectutils import Epoch, Iteration
from biobb_md.gromacs.pdb2gmx import Pdb2gmx
from biobb_md.gromacs.solvate import Solvate
from molecules import (
    GROComplex,
    ZipTopology,
    GROStructure,
    PDBStructure,
)
from typing import Dict
from primitives import launch_biobb

# DEPRECATED
def minimize_iterations(work_dir: WorkProject, iter_name: str, maxsol: int = 1):
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
    pdb_to_gro_zip = Pdb2gmx(
        input_pdb_path=str(this_iter.complex.pdb.file.path),
        output_gro_path=dry_gro,
        output_top_zip_path=dry_top_zip,
        properties=props,
    )
    launch_biobb(pdb_to_gro_zip)

    # SOLVATE
    wet_gro = str(current_dir / ("wet_" + system_name + ".gro"))
    wet_top_zip = str(current_dir / "wet_topol.zip")
    props = {"gmx_path": str(config["md"]["gmx_bin"]), "dev": f"-maxsol {maxsol}"}
    solvatador = Solvate(
        input_solute_gro_path=dry_gro,
        output_gro_path=wet_gro,
        input_top_zip_path=dry_top_zip,
        output_top_zip_path=wet_top_zip,
        properties=props,
    )
    launch_biobb(solvatador)

    a = ZipTopology.from_path_with_chains(
        wet_top_zip,
        target_chains=work_dir.config["target"]["chainID"],
        binder_chains=work_dir.config["binder"]["chainID"],
    )
    b = GROStructure.from_path(wet_gro)
    c = PDBStructure.from_path(this_iter.complex.pdb.file.path)
    this_iter.complex = GROComplex(system_name, current_dir, c, a, b)
    this_iter.complex.update_all_ndxs()
