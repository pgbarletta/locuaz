import warnings
from collections.abc import Iterable
from pathlib import Path
from typing import Tuple

import MDAnalysis as mda
import numpy as np
from biobb_analysis.gromacs.gmx_image import GMXImage

from locuaz.complex import GROComplex
from locuaz.fileutils import FileHandle
from locuaz.molecules import PDBStructure, read_ndx
from locuaz.primitives import launch_biobb


def fix_box(uni: mda.Universe, *, target_indices: Iterable, binder_indices: Iterable,
            non_protein_indices: Iterable) -> Tuple[bool, mda.Universe]:
    """
    fix_box() wraps the atoms of a system around its box and then centers them.
    From: Baptista, A. M.; da Rocha, L.; Campos, S. R. R. FixBox: A General
    Algorithm to Fix Molecular Systems in Periodic Boxes. J. Chem. Inf. Model.
    2022, 62 (18), 4435–4447. https://doi.org/10.1021/acs.jcim.2c00823.

    Args:
        uni (mda.Universe):
        target_indices (Iterable):
        binder_indices (Iterable):
        non_protein_indices (Iterable):

    Returns:
        Tuple[bool, mda.Universe]: boolean indicating whether there are target or
         binder atoms outside the box and the universe with the wrapped coords.
    """
    H = get_matrix(uni.dimensions)
    inv_H = np.linalg.inv(H)
    centro = np.sum(H / 2, axis=0)

    # s_positions = (uni.atoms.positions * 0.1 - centro) @ inv_H  # type:ignore
    s_positions = (uni.atoms.positions - centro) @ inv_H  # type:ignore

    # Reassemble complex
    # TODO: Add chain indices for target and binder, so this works with multi-chain.
    # If not, binder_indices and target_indices will have atoms from different chains
    # that may be on different images. Then we can assume all chains from target and binder have to
    # be clustered together separately first, and then we cluster the 2 groups.
    min_distances = []
    lista_idx_mini = []
    for i, k in enumerate(binder_indices):
        ds_i = s_positions[target_indices] - s_positions[binder_indices][i, :]
        ds_i_imaged = ds_i - np.floor(ds_i + 0.5)
        dist_i_imaged = np.sum((ds_i_imaged @ H) ** 2, axis=1)
        idx_min = np.argmin(dist_i_imaged)
        mini = dist_i_imaged[idx_min]

        lista_idx_mini.append(idx_min)
        min_distances.append(mini)

    binder_closest = np.argmin(min_distances)
    target_closest = lista_idx_mini[binder_closest]

    ds_i = (
            s_positions[target_indices][target_closest]
            - s_positions[binder_indices][binder_closest, :]
    )
    box_displacement = np.floor(ds_i + 0.5)
    s_positions[binder_indices] = s_positions[binder_indices] + box_displacement

    # Center complex
    protein_coords = s_positions[np.append(target_indices, binder_indices)]
    box_x = (np.min(protein_coords[:, 0]) + np.max(protein_coords[:, 0])) / 2
    box_y = (np.min(protein_coords[:, 1]) + np.max(protein_coords[:, 1])) / 2
    box_z = (np.min(protein_coords[:, 2]) + np.max(protein_coords[:, 2])) / 2
    box = [box_x, box_y, box_z]
    s_positions = s_positions - box

    # Rewrap solvent in box
    wat_oxygens = [atm for atm in uni.atoms[non_protein_indices] if
                   atm.element == "O" and atm.resname in {"SOL", "WAT"}]  # type: ignore
    waters = [atm.residue for atm in wat_oxygens]  # type: ignore

    assert len(waters) == len(wat_oxygens), "This should not happen fix_box() failed. "\
        f"Number of O water atoms: {len(wat_oxygens)}, number of water molecules: {len(waters)}"

    for wat, oxy in zip(waters, wat_oxygens):
        wat_atm_indices = wat.atoms.indices
        O_xyz = s_positions[oxy.index]
        wrapped_O_xyz = np.floor(O_xyz + 0.5)
        for i in wat_atm_indices:
            s_positions[i] -= wrapped_O_xyz

    # Re-wrap non-protein that aren't solvent. This should be just ions.
    ions_residues = {
        atm.residue for atm in uni.atoms[non_protein_indices] if atm.resname not in {"SOL", "WAT"}  # type: ignore
    }
    for ion in ions_residues:
        ion_atm_indices = ion.atoms.indices
        # Wrap them around using the first atom of the residue.
        ion_xyz = s_positions[ion_atm_indices[0]]
        wrapped_ion_xyz = np.floor(ion_xyz + 0.5)
        for i in ion_atm_indices:
            s_positions[i] -= wrapped_ion_xyz

    uni.atoms.positions = ((s_positions @ H) + centro)  # type: ignore

    # Check if successful
    n_outside_box_target = np.sum(np.floor(s_positions[target_indices] + 0.5))
    n_outside_box_binder = np.sum(np.floor(s_positions[binder_indices] + 0.5))
    all_in = (n_outside_box_target + n_outside_box_binder) == 0

    return all_in, uni


def fix_box_cpx(
    cpx: GROComplex, out_path: Path, gmx_bin: str = "gmx"
) -> Tuple[bool, PDBStructure]:
    """
    fix_box_cpx(): fix_box() wrapper for the protocol.
    Args:
        cpx (GROComplex): GROComplex
        out_path (Path): output PDB
        gmx_bin (str): GROMACS binary executable name

    Returns:
        Tuple[bool, PDBStructure]: boolean indicating whether there are target or
         binder atoms outside the box and PDBStructure of the input `out_path`.
    """

    # First, make sure the protein is whole. This is a requirement for FixBox to work
    out_path = Path(out_path)
    whole_pdb = Path(out_path.parent, "whole.pdb")
    make_whole = GMXImage(
        input_traj_path=str(cpx.gro),
        input_top_path=str(cpx.tpr),
        output_traj_path=str(whole_pdb),
        properties={
            "binary_path": gmx_bin,
            "fit_selection": "Protein",
            "center_selection": "System",
            "output_selection": "System",
            "pbc": "whole",
            "center": False,
        },
    )
    launch_biobb(make_whole)

    # Now, run FixBox
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        u = mda.Universe(str(whole_pdb), in_memory=True)
    indices = read_ndx(cpx.ndx.path)
    all_in, u = fix_box(u, target_indices=indices["target"],
                        binder_indices=indices["binder"],
                        non_protein_indices=indices["Non-Protein"])

    # Write out the modified Universe
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        u.atoms.write(out_path)  # type: ignore

    # Remove temporary files
    whole_pdb.unlink()

    return all_in, PDBStructure(FileHandle(out_path))


def get_matrix(dimensions):
    x, y, z, a, b, c = dimensions
    # x /= 10
    # y /= 10
    # z /= 10
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
