import shutil as sh
from collections.abc import Iterable
from pathlib import Path
from typing import Dict, Tuple, Optional, Union

import MDAnalysis as mda
import numpy as np
from biobb_analysis.gromacs.gmx_trjconv_str import GMXTrjConvStr
from biobb_gromacs.gromacs.editconf import Editconf
from biobb_gromacs.gromacs.genion import Genion
from biobb_gromacs.gromacs.grompp import Grompp
from biobb_gromacs.gromacs.pdb2gmx import Pdb2gmx

from amberutils import create_tleap_script, fix_pdb, run_tleap, amb_to_gmx
from fileutils import FileHandle, copy_to
from molecules import PDBStructure, GROStructure, ZipTopology
from primitives import launch_biobb


def get_gro_ziptop_from_pdb(
        *,
        pdb: PDBStructure,
        target_chains: Iterable,
        binder_chains: Iterable,
        md_config: Dict,
        add_ions: bool = False,
) -> Tuple[PDBStructure, GROStructure, ZipTopology]:
    """get_gro_ziptop_from_pdb does a pdb2gmx from the PDB and tries to keep
    the system neutral, which may alter the topology so a new PDB will be
    written with the same name as the original, which will be backed up by GROMACS.

    Args:
        pdb (PDBStructure): input PDB
        target_chains (Iterable): these will be used to construct the ZipTopology
        binder_chains (Iterable): these will be used to construct the ZipTopology

    Returns:
        Tuple[PDBStructure, GROStructure, ZipTopology]: Proper, nice, system.
    """
    gmx_bin: str = md_config.get("gmx_bin", "gmx")
    water_type: str = md_config.get("water_type", "tip3p")
    force_field: str = md_config.get("force_field", "amber99sb-ildn")

    local_dir = pdb.file.path.parent
    name = pdb.name

    # Generate the first set of GROStructure and ZipTopology
    props = {
        "binary_path": gmx_bin,
        "water_type": water_type,
        "force_field": force_field,
        "ignh": True,
    }
    pre_gro_fn = local_dir / ("pre" + name + ".gro")
    pre_top_fn = local_dir / "pre_topol.zip"
    pdb_to_gro_zip = Pdb2gmx(
        input_pdb_path=str(pdb.file),
        output_gro_path=str(pre_gro_fn),
        output_top_zip_path=str(pre_top_fn),
        properties=props,
    )
    launch_biobb(pdb_to_gro_zip)

    if md_config.get("use_box", False):
        # Set the box dimensions
        box: Optional[str] = md_config.get("box")
        if box:
            props = {
                "distance_to_molecule": None,
                "box_type": "triclinic",
                "dev": f"-box {box}",
            }
        else:
            dist_to_box: float = md_config.get("dist_to_box", 1.0)
            box_type: str = md_config.get("box_type", "triclinic")
            props = {
                "distance_to_molecule": dist_to_box,
                "box_type": box_type,
            }

        box_gro_fn = local_dir / ("box" + name + ".gro")
        set_box = Editconf(
            input_gro_path=str(pre_gro_fn),
            output_gro_path=str(box_gro_fn),
            properties=props,
        )
        launch_biobb(set_box)

        pre_gro_fn = copy_to(FileHandle(box_gro_fn), local_dir, f"pre{name}.gro").path
        box_gro_fn.unlink()

    if add_ions:
        # Build a .tpr for genion
        gen_tpr_fn = local_dir / ("genion_" + name + ".tpr")
        grompepe = Grompp(
            input_gro_path=str(pre_gro_fn),
            input_top_zip_path=str(pre_top_fn),
            output_tpr_path=str(gen_tpr_fn),
            properties={"binary_path": gmx_bin, "maxwarn": 2},
        )
        launch_biobb(grompepe)

        # Add ions
        gro_fn = local_dir / (name + ".gro")
        top_fn = local_dir / (name + ".zip")
        genio = Genion(
            input_tpr_path=str(gen_tpr_fn),
            input_top_zip_path=str(pre_top_fn),
            output_gro_path=str(gro_fn),
            output_top_zip_path=str(top_fn),
            properties={"binary_path": gmx_bin, "neutral": True, "concentration": 0.0},
        )
        launch_biobb(genio)
        # Remove temporary file
        # gen_tpr_fn.unlink()
    else:
        gro_fn = copy_to(FileHandle(pre_gro_fn), local_dir, name + ".gro").path
        top_fn = copy_to(FileHandle(pre_top_fn), local_dir, name + ".zip").path

    # Build a temporary tpr file for the next step
    temp_tpr_fn = local_dir / ("temp_" + name + ".tpr")
    grompepe = Grompp(
        input_gro_path=str(gro_fn),
        input_top_zip_path=str(top_fn),
        output_tpr_path=str(temp_tpr_fn),
        properties={
            "binary_path": gmx_bin,
            "maxwarn": 2,
        },
    )
    launch_biobb(grompepe)

    # Get PDB from GRO file. Gromacs should back up the older input PDB, maybe?
    pdb_fn = local_dir / (name + ".pdb")
    trjconv = GMXTrjConvStr(
        input_structure_path=str(gro_fn),
        input_top_path=str(temp_tpr_fn),
        output_str_path=str(pdb_fn),
        properties={"binary_path": gmx_bin},
    )
    launch_biobb(trjconv)

    # Remove temporary files
    temp_tpr_fn.unlink()
    pre_gro_fn.unlink()
    pre_top_fn.unlink()

    zip_top = ZipTopology.from_path_with_chains(
        top_fn, target_chains=target_chains, binder_chains=binder_chains
    )

    return PDBStructure.from_path(pdb_fn), GROStructure.from_path(gro_fn), zip_top


def get_gro_ziptop_from_pdb_tleap(
    *,
    pdb: PDBStructure,
    target_chains: Iterable,
    binder_chains: Iterable,
) -> Tuple[PDBStructure, GROStructure, ZipTopology]:
    """get_gro_ziptop_from_pdb_tleap runs .

    Args:
        pdb (PDBStructure): input PDB
        target_chains (Iterable): these will be used to construct the ZipTopology
        binder_chains (Iterable): these will be used to construct the ZipTopology

    Returns:
        Tuple[PDBStructure, GROStructure, ZipTopology]: Proper, nice, system.
    """
    pdb_path = Path(pdb)
    local_dir = pdb_path.parent
    name = pdb.name

    tleap_script = create_tleap_script(local_dir, name)

    # Backup the PDB before runing pdb4amber
    pre_fix_pdb = local_dir / f"preAmberPDBFixer_{pdb_path.stem}.pdb"
    sh.move(pdb_path, pre_fix_pdb)

    pdb = fix_pdb(pre_fix_pdb, pdb_path)

    # Backup the PDB before runing tleap
    pre_tleap_pdb = local_dir / f"pretleap_{pdb_path.stem}.pdb"
    sh.copy(Path(pdb), pre_tleap_pdb)
    pdb, prmtop, rst = run_tleap(tleap_script, name)

    gro, ziptop = amb_to_gmx(
        name,
        pdb,
        prmtop,
        rst,
        target_chains=target_chains,
        binder_chains=binder_chains,
    )

    return pdb, gro, ziptop


def fix_wat_naming(
    pdb_in: Union[PDBStructure, Path], pdb_out: Path, *, use_tleap: bool = False
) -> PDBStructure:

    u = mda.Universe(str(pdb_in))
    if use_tleap:
        sol_atomgroup = u.select_atoms("resname SOL")
        if len(sol_atomgroup) > 0:
            for r in sol_atomgroup.residues:
                r.resname = "WAT"
                r.atoms.names = np.array(["O", "H1", "H2"])
    else:
        wat_atomgroup = u.select_atoms("resname WAT")
        if len(wat_atomgroup) > 0:
            for r in wat_atomgroup.residues:
                r.resname = "SOL"
                r.atoms.names = np.array(["OW", "HW1", "HW2"])

    u.atoms.write(str(pdb_out))  # type: ignore

    return PDBStructure.from_path(pdb_out)
