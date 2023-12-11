import shutil as sh
import warnings
from collections.abc import Iterable
from pathlib import Path
from typing import Dict, Tuple, Union

import MDAnalysis as mda
import numpy as np
from biobb_analysis.gromacs.gmx_trjconv_str import GMXTrjConvStr
from biobb_gromacs.gromacs.genion import Genion
from biobb_gromacs.gromacs.pdb2gmx import Pdb2gmx

from locuaz.amberutils import create_tleap_script, fix_pdb, run_tleap, amb_to_gmx
from locuaz.molecules import PDBStructure, GROStructure, ZipTopology, get_tpr
from locuaz.primitives import launch_biobb

from multiprocessing import Lock
lock = Lock()

def get_gro_ziptop_from_pdb(
        *,
        pdb: PDBStructure,
        target_chains: Iterable,
        binder_chains: Iterable,
        md_config: Dict,
        add_ions: bool = False,
) -> Tuple[PDBStructure, GROStructure, ZipTopology]:
    """
    Does a pdb2gmx from the PDB and tries to keep the system neutral,
    which may alter the topology so a new PDB will be written with the same name as the original,
    which will be backed up by GROMACS. No 'editconf' or 'solvate' is run.
    The system is assumed to have solvent and a CRYST1 record with box info.

    Parameters
    ----------
    pdb : PDBStructure
        input PDB
    target_chains : Iterable
        target chainIDs. Will be used to construct the ZipTopology
    binder_chains : Iterable
        binder chainIDs. Will be used to construct the ZipTopology
    md_config : Dict
        user's config input for MD
    add_ions
        whether to neutralize with genion or not. Useful after a mutation was
        done on the input PDB.
    Returns
    -------
    pdb, gro, ziptop : Tuple[PDBStructure, GROStructure, ZipTopology]
        Proper, nice, system.
    """
    gmx_bin: str = "gmx"
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
    pre_ion_gro_fn = local_dir / f"pre_ion_{name}.gro"
    pre_ion_top_fn = local_dir / f"pre_ion_{name}.zip"
    pdb_to_gro_zip = Pdb2gmx(
        input_pdb_path=str(pdb.file),
        output_gro_path=str(pre_ion_gro_fn),
        output_top_zip_path=str(pre_ion_top_fn),
        properties=props,
    )
    with lock:
        launch_biobb(pdb_to_gro_zip)

    gro_fn = local_dir / f"{name}.gro"
    top_fn = local_dir / f"{name}.zip"
    if add_ions:
        # Build a temporary .tpr for genion
        temp_tpr_fn = get_tpr(gro=pre_ion_gro_fn, top=pre_ion_top_fn)

        # Re-add ions as necessary. SOL group will be continuous, so gmx genion won't complain.
        genio = Genion(
            input_tpr_path=str(temp_tpr_fn),
            input_top_zip_path=str(pre_ion_top_fn),
            output_gro_path=str(gro_fn),
            output_top_zip_path=str(top_fn),
            properties={
                "binary_path": gmx_bin,
                "neutral": True,
                "concentration": 0.0,
            },
        )
        launch_biobb(genio)
    else:
        sh.move(pre_ion_gro_fn, gro_fn)
        sh.move(pre_ion_top_fn, top_fn)

    # Build a temporary tpr file for the next step
    temp_tpr_fn = get_tpr(gro=gro_fn, top=top_fn)

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
    """

    Parameters
    ----------
    pdb : PDBStructure
        input PDB
    target_chains : Iterable
        target chainIDs. Will be used to construct the ZipTopology
    binder_chains : Iterable
        binder chainIDs. Will be used to construct the ZipTopology

    Returns
    -------
    pdb, gro, ziptop : Tuple[PDBStructure, GROStructure, ZipTopology]
        Proper, nice, system.
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
    # Will get some warnings due to new waters added by gmx solvate not having the element column
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
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

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        u.atoms.write(str(pdb_out))  # type: ignore

    return PDBStructure.from_path(pdb_out)


def set_posres(in_posres: Path, out_posres: Path, restraint: float) -> None:
    with open(in_posres, "r") as file:
        txt_in_posres = file.readlines()

    with open(out_posres, "w") as file:
        for linea in txt_in_posres:
            new_line = linea
            if linea[0] not in {';', '[', '\n'}:
                # This line should look something like: '    13     1  1000  1000  1000'
                try:
                    atm = int(linea.split()[0])
                    new_line = f"{atm:6d}{1:6d}{restraint:6d}{restraint:6d}{restraint:=6d}\n"
                except ValueError:
                    # well, apparently it didn't
                    pass
            file.write(new_line)
    return
