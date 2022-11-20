from pathlib import Path
import numpy as np
from typing import Any, Tuple
from functools import singledispatch

import MDAnalysis as mda
from biobb_analysis.gromacs.gmx_image import GMXImage
from biobb_analysis.gromacs.gmx_trjconv_str import GMXTrjConvStr

from fileutils import FileHandle
from primitives import launch_biobb
from molecules import PDBStructure, read_ndx
from complex import AbstractComplex, GROComplex


@singledispatch
def split_solute_and_solvent(complex: AbstractComplex, gmx_bin: str) -> Any:
    raise NotImplementedError


@split_solute_and_solvent.register
def _(complex: GROComplex, gmx_bin: str) -> Tuple[PDBStructure, PDBStructure]:
    """prepare_old_iter extract 2 PDBs from an input pdb, one with the protein
    and the other with the water and ions.

    Args:
        complex (Complex): a complex object with a PDB and a TPR file.
    """

    # Protein
    nonwat_pdb_fn = Path(complex.dir) / ("nonwat_" + complex.name + ".pdb")
    get_protein = GMXTrjConvStr(
        input_structure_path=str(complex.pdb.file.path),
        input_top_path=str(complex.tpr.file.path),
        input_index_path=str(complex.ndx.path),
        output_str_path=str(nonwat_pdb_fn),
        properties={"gmx_path": gmx_bin, "selection": "Protein"},
    )
    launch_biobb(get_protein)
    nonwat_pdb = PDBStructure.from_path(nonwat_pdb_fn)

    # Water and ions
    wation_pdb_fn = Path(complex.dir) / ("wation_" + complex.name + ".pdb")
    get_water_ions = GMXTrjConvStr(
        input_structure_path=str(complex.pdb.file.path),
        input_top_path=str(complex.tpr.file.path),
        input_index_path=str(complex.ndx.path),
        output_str_path=str(wation_pdb_fn),
        properties={"gmx_path": gmx_bin, "selection": "Non-Protein"},
    )
    launch_biobb(get_water_ions)

    wation_pdb = PDBStructure.from_path(wation_pdb_fn)

    return nonwat_pdb, wation_pdb


def get_matrix(dimensions):
    x, y, z, a, b, c = dimensions
    x /= 10
    y /= 10
    z /= 10
    H = np.zeros((3, 3))
    H[0, 0] = x
    if a == 90.0 and b == 90.0 and c == 90.0:
        H[1, 1] = y
        H[2, 2] = z
    else:
        a = np.deg2rad(a)
        b = np.deg2rad(b)
        c = np.deg2rad(c)
        H[1][0] = y * np.cos(c)
        H[1][1] = y * np.sin(c)
        H[2][0] = z * np.cos(b)
        H[2][1] = z * (np.cos(a) - np.cos(b) * np.cos(c)) / np.sin(c)
        H[2][2] = np.sqrt(z * z - H[2][0] ** 2 - H[2][1] ** 2)
    return H


def fix_box(cpx: GROComplex, out_path: Path, gmx_bin: str = "gmx") -> PDBStructure:

    whole_pdb = Path(out_path.parent, "whole.pdb")
    make_whole = GMXImage(
        input_traj_path=str(cpx.gro),
        input_top_path=str(cpx.tpr),
        output_traj_path=str(whole_pdb),
        properties={
            "gmx_path": gmx_bin,
            "fit_selection": "Protein",
            "center_selection": "System",
            "output_selection": "System",
            # "ur": "compact",
            "pbc": "whole",
            "center": False,
        },
    )
    launch_biobb(make_whole)

    u = mda.Universe(whole_pdb)
    H = get_matrix(u.dimensions)
    inv_H = np.linalg.inv(H)
    centro = np.sum(H / 2, axis=1)

    s_positions = (u.atoms.positions * 0.1 - centro) @ inv_H  # type:ignore
    indices = read_ndx(cpx.ndx.path)

    # Reassemble complex
    min_distances = []
    lista_idx_mini = []
    for i, k in enumerate(indices["binder"]):
        ds_i = s_positions[indices["target"]] - s_positions[indices["binder"]][i, :]
        ds_i_imaged = ds_i - np.floor(ds_i + 0.5)
        dist_i_imaged = np.sum((ds_i_imaged @ H) ** 2, axis=1)
        idx_min = np.argmin(dist_i_imaged)
        mini = dist_i_imaged[idx_min]

        lista_idx_mini.append(idx_min)
        min_distances.append(mini)

    binder_closest = np.argmin(min_distances)
    target_closest = lista_idx_mini[binder_closest]

    ds_i = (
        s_positions[indices["target"]][target_closest]
        - s_positions[indices["binder"]][binder_closest, :]
    )
    box_displacement = np.floor(ds_i + 0.5)
    s_positions[indices["binder"]] = s_positions[indices["binder"]] + box_displacement

    # Center complex
    protein_coords = s_positions[np.append(indices["target"], indices["binder"])]
    box_x = (np.min(protein_coords[:, 0]) + np.max(protein_coords[:, 0])) / 2
    box_y = (np.min(protein_coords[:, 1]) + np.max(protein_coords[:, 1])) / 2
    box_z = (np.min(protein_coords[:, 2]) + np.max(protein_coords[:, 2])) / 2
    box = [box_x, box_y, box_z]
    s_positions = s_positions - box

    # Rewrap solvent in box
    waters = {atm.residue for atm in u.atoms[indices["Non-Protein"]] if atm.resname == "SOL"}  # type: ignore
    wat_oxygens = {atm for atm in u.atoms[indices["Non-Protein"]] if atm.resname == "SOL" and atm.element == "O"}  # type: ignore

    for wat, oxy in zip(waters, wat_oxygens):
        wat_atm_indices = wat.atoms.indices
        O_xyz = s_positions[oxy.index]
        wrapped_O_xyz = np.floor(O_xyz + 0.5)
        for i in wat_atm_indices:
            s_positions[i] -= wrapped_O_xyz

    # Rewrap non-protein that aren't solvent. This should be just ions.
    ions_residues = {
        atm.residue for atm in u.atoms[indices["Non-Protein"]] if atm.resname != "SOL"  # type: ignore
    }
    for ion in ions_residues:
        ion_atm_indices = ion.atoms.indices
        # Wrap them around using the first atom of the residue.
        ion_xyz = s_positions[ion_atm_indices[0]]
        wrapped_ion_xyz = np.floor(ion_xyz + 0.5)
        for i in ion_atm_indices:
            s_positions[i] -= wrapped_ion_xyz

    u.atoms.positions = ((s_positions @ H) + centro) * 10  # type: ignore
    u.atoms.write(out_path)  # type: ignore

    # Remove temporary files
    whole_pdb.unlink()

    return PDBStructure(FileHandle(out_path))
