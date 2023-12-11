from abc import ABCMeta
from collections.abc import Iterable
from functools import singledispatch
from itertools import chain
from pathlib import Path
from typing import List, Dict, Tuple, Any, Union

import MDAnalysis as mda
import numpy as np
from attrs import define, field
from biobb_analysis.gromacs.gmx_trjconv_str import GMXTrjConvStr
from biobb_gromacs.gromacs.gmxselect import Gmxselect
from biobb_gromacs.gromacs.grompp import Grompp
from biobb_gromacs.gromacs.pdb2gmx import Pdb2gmx

from locuaz.fileutils import FileHandle, DirHandle, copy_to, update_header
from locuaz.primitives import launch_biobb


@define
class AbstractFileObject(metaclass=ABCMeta):
    file: FileHandle = field(converter=FileHandle)  # type: ignore
    name: str = field(init=False)
    ext: str = field(init=False)

    def __attrs_post_init__(self):
        self.name, self.ext = self.file.path.name.split(".")

    def unlink(self) -> None:
        self.file.unlink()

    def __str__(self) -> str:
        return str(self.file)

    def __fspath__(self) -> str:
        return self.__str__()

    @classmethod
    def from_path(cls, path: Any):
        file = FileHandle(path)
        return cls(file)


@define
class TPRFile(AbstractFileObject):
    pass


@define
class Structure(AbstractFileObject):
    pass


@define
class Trajectory(AbstractFileObject):
    pass


@define
class PDBStructure(Structure):
    """PDBStructure

    Args:
        file (FileHandle): path to the file
    """

    # TODO: turn __init__ to __attrs_post_init__
    def __init__(self, file: FileHandle):
        super().__init__(file)
        if self.ext != "pdb":
            raise ValueError("Input file's extension for PDBStructure should be .pdb")

    def get_cryst1_record(self) -> str:
        with open(self.file.path, "r") as file:
            for linea in file:
                if linea[0:6] == "CRYST1":
                    cryst1_record = linea
                    break
            else:
                raise ValueError(f"{self} has no CRYST1 record.")
        return cryst1_record

    def set_cryst1_record(self, cryst1_record: str) -> None:
        texto = []
        with open(self.file.path, "r") as file:
            for linea in file:
                if linea[0:6] == "CRYST1":
                    continue
                texto.append(linea)
        with open(self.file.path, "w") as file:
            if cryst1_record[-1] == "\n":
                file.write(cryst1_record)
            else:
                file.write(cryst1_record + "\n")
            [file.write(linea) for linea in texto]


@define
class GROStructure(Structure):
    """GROStructure

    Args:
        file (FileHandle): path to the file
    """

    # TODO: turn into __attrs_post_init__
    def __init__(self, file: FileHandle):
        super().__init__(file)
        if self.ext != "gro":
            raise ValueError("Input file's extension for PDBStructure should be .gro")


@define
class Topology(AbstractFileObject):
    target_chains: Tuple[str] = field(init=False)
    binder_chains: Tuple[str] = field(init=False)

    selection_target: str = field(init=False)
    selection_binder: str = field(init=False)
    selection_complex: str = field(init=False)
    selection_not_complex: str = field(init=False)
    selection_water: str = field(init=False)

    @classmethod
    def from_path_with_chains(
        cls, path: Path, *, target_chains: Iterable, binder_chains: Iterable
    ) -> "Topology":
        self = super().from_path(path)
        self.target_chains = tuple(target_chains)
        self.binder_chains = tuple(binder_chains)

        self.selection_target = " or ".join(
            [f"segid {chainID}" for chainID in target_chains]
        )
        self.selection_binder = " or ".join(
            [f"segid {chainID}" for chainID in binder_chains]
        )
        self.selection_complex = " or ".join(
            ["protein", self.selection_target, self.selection_binder]
        )
        self.selection_not_complex = " and not ".join(
            ["not protein", self.selection_target, self.selection_binder]
        )
        self.selection_water = "resname WAT or resname SOL"

        return self


@define
class AmberTopology(Topology):
    pass


@define
class ZipTopology(Topology):
    @classmethod
    def from_path_with_chains(
        cls, path: Path, *, target_chains: Iterable, binder_chains: Iterable
    ) -> "ZipTopology":
        self = super().from_path(path)
        self.target_chains = tuple(target_chains)
        self.binder_chains = tuple(binder_chains)

        self.selection_target = " or ".join(
            [f"segid {chainID}" for chainID in target_chains]
        )
        self.selection_binder = " or ".join(
            [f"segid {chainID}" for chainID in binder_chains]
        )
        self.selection_complex = " or ".join(
            ["protein", f"({self.selection_target})", f"({self.selection_binder})"]
        )
        self.selection_not_complex = " and not ".join(
            ["not protein", f"({self.selection_target})", f"({self.selection_binder})"]
        )
        self.selection_water = "resname WAT or resname SOL"

        return self


@define
class TrrTrajectory(Trajectory):
    pass


@define
class XtcTrajectory(Trajectory):
    pass


def get_tpr(
    *,
    gro: Union[Path, GROStructure],
    top: Union[Path, ZipTopology],
    gmx_bin: str = "gmx",
) -> TPRFile:
    """
    get_tpr(): build a quick tpr, useful for simple runs.
    Args:
        gro:
        top:
        gmx_bin:

    Returns:
        TPRFile
    """
    gro_fn = Path(gro)
    tpr_fn = gro_fn.parent / (gro_fn.stem + ".tpr")
    grompepe = Grompp(
        input_gro_path=str(gro_fn),
        input_top_zip_path=str(top),
        output_tpr_path=str(tpr_fn),
        properties={"binary_path": gmx_bin, "simulation_type": "index"},
    )
    launch_biobb(grompepe)

    return TPRFile.from_path(tpr_fn)


def get_pdb_tpr(
    *,
    gro: GROStructure,
    top: ZipTopology,
    gmx_bin: str = "gmx",
) -> Tuple[PDBStructure, TPRFile]:

    tpr = get_tpr(gro=gro, top=top, gmx_bin=gmx_bin)

    pdb_fn = Path(gro.file.path.parent) / (gro.name + ".pdb")
    trjconv = GMXTrjConvStr(
        input_structure_path=str(gro.file.path),
        input_top_path=str(tpr.file.path),
        output_str_path=str(pdb_fn),
        properties={"binary_path": gmx_bin},
    )
    launch_biobb(trjconv)

    return PDBStructure.from_path(pdb_fn), tpr


def fix_pdb_gro(
    pdb_in: PDBStructure, pdb_out_path: Path, gmx_bin: str = "gmx"
) -> PDBStructure:
    # This is quite cumbersome.
    wrk_dir = pdb_in.file.path.parent

    # First, use pdb2gmx to fix the PDB file. There's no way of ignoring
    # Hydrogens and getting a good output PDB, so the output has to be
    # a .gro file
    props = {
        "binary_path": gmx_bin,
        "water_type": "tip3p",
        "force_field": "amber99sb-ildn",
        "ignh": True,
    }

    temp_gro = wrk_dir / "temp.gro"
    temp_zip = wrk_dir / "temp.zip"
    pdb_fix = Pdb2gmx(
        input_pdb_path=str(pdb_in.file.path),
        output_gro_path=str(temp_gro),
        output_top_zip_path=str(temp_zip),
        properties=props,
    )
    launch_biobb(pdb_fix)

    # Now, get a tpr file from the new GRO and zip files
    temp_tpr = wrk_dir / "temp.tpr"
    grompepe = Grompp(
        input_gro_path=str(temp_gro),
        input_top_zip_path=str(temp_zip),
        output_tpr_path=str(temp_tpr),
        properties={"binary_path": gmx_bin, "maxwarn": 1},
    )
    launch_biobb(grompepe)

    # Finally, get a PDB file from the GRO file
    trjconv = GMXTrjConvStr(
        input_structure_path=str(temp_gro),
        input_top_path=str(temp_tpr),
        output_str_path=str(pdb_out_path),
        properties={"binary_path": gmx_bin},
    )
    launch_biobb(trjconv)

    temp_gro.unlink()
    temp_zip.unlink()
    temp_tpr.unlink()

    return PDBStructure.from_path(pdb_out_path)


def catenate_pdbs(
    *pdbs: PDBStructure, pdb_out_path: Path, gmx_bin: str = "gmx"
) -> PDBStructure:
    lineas = []
    for pdb in pdbs:
        with open(pdb.file.path, "r") as file:
            for linea in file:
                if (
                    linea[0:4] == "ATOM"
                    or linea[0:3] == "TER"
                    or linea[0:6] == "HETATM"
                ):
                    lineas.append(linea)

    texto = chain.from_iterable(lineas)
    temp_pdb = pdb_out_path.parent / "temp.pdb"
    with open(temp_pdb, "w") as file:
        [file.write(linea) for linea in texto]
        file.write("END")

    # Now, use gromacs to fix the PDB.
    out_pdb = fix_pdb_gro(PDBStructure.from_path(temp_pdb), pdb_out_path, gmx_bin)
    temp_pdb.unlink()

    return out_pdb


def write_chain_ndx(
    *, pdb: PDBStructure, chains: Iterable, selname: str, gmx_bin: str = "gmx"
) -> FileHandle:
    selection = " or ".join([f"chain {chainID}" for chainID in chains])

    wrk_dir = pdb.file.path.parent
    ndx_fn = Path(wrk_dir) / f"{selname}.ndx"
    selector = Gmxselect(
        input_structure_path=str(pdb),
        output_ndx_path=str(ndx_fn),
        properties={"selection": selection},
    )
    launch_biobb(selector)

    ndx = FileHandle(ndx_fn)
    update_header(ndx, f"[ {selname} ]\n")

    return ndx


def generate_ndx(
    name: str,
    *,
    pdb: PDBStructure,
    top: Topology,
) -> FileHandle:
    try:
        uni_pdb = mda.Universe(str(pdb))
        ndx_file = Path(Path(pdb).parent, f"{name}.ndx")

        # sel_target = " or ".join([f"segid {chainID}" for chainID in target_chains])
        uni_pdb.select_atoms(top.selection_target).write(
            ndx_file, name="target", mode="w"
        )

        # sel_binder = " or ".join([f"segid {chainID}" for chainID in binder_chains])
        uni_pdb.select_atoms(top.selection_binder).write(
            ndx_file, name="binder", mode="a"
        )

        # TODO: this won't work with glyco mods and stuff
        # sel_protein = " or ".join(["protein", sel_target, sel_binder])
        uni_pdb.select_atoms(top.selection_complex).write(
            ndx_file, name="Protein", mode="a"
        )

        # sel_not_target = " and ".join([f"segid {chainID}" for chainID in target_chains])
        # sel_not_binder = " and ".join([f"segid {chainID}" for chainID in binder_chains])
        # sel_not_protein = " and not ".join(
        #     ["not protein", sel_not_target, sel_not_binder]
        # )
        uni_pdb.select_atoms(top.selection_not_complex).write(
            ndx_file, name="Non-Protein", mode="a"
        )
        uni_pdb.select_atoms("all").write(ndx_file, name="sistema", mode="a")

        return FileHandle(ndx_file)

    except Exception as e:
        raise e


def read_ndx(ndx_path: Path) -> Dict:
    """
    read_ndx(): reads GROMACS .ndx file. Watch out, it also performs 0-index
    correction. That is, it substracts 1 from the atom indices.
    Args:
        ndx_path(Path):

    Returns:
        Dict: 0-indexed selections from the .ndx file
    """
    indices: List[str] = []
    group_name = None
    grupos: Dict[str, Any] = {}
    with open(ndx_path, "r") as ndx_file:
        for line in ndx_file:
            if line.startswith("["):
                if group_name:
                    grupos[group_name] = np.array(indices, dtype=int) - 1
                group_name = line[1:-2].strip()
                indices = []
            else:
                indices += line.split()
    return grupos


@singledispatch
def copy_mol_to(obj, dir_path: Union[Path, DirHandle], name=None):
    raise NotImplementedError


@copy_mol_to.register
def _(obj: Structure, dir_path: Union[Path, DirHandle], name=None):
    new_file = copy_to(obj.file, dir_path, name=name)
    return Structure(new_file)


@copy_mol_to.register
def _(obj: ZipTopology, dir_path: Union[Path, DirHandle], name=None):
    new_file = copy_to(obj.file, dir_path, name=name)
    ZipTopology.from_path_with_chains(
        new_file.path, target_chains=obj.target_chains, binder_chains=obj.binder_chains
    )


@copy_mol_to.register
def _(obj: Trajectory, dir_path: Union[Path, DirHandle], name=None):
    new_file = copy_to(obj.file, dir_path, name=name)
    return Trajectory(new_file)


def try_copy_to(obj, dir_path: Union[Path, DirHandle], name=None) -> Any:
    try:
        new_obj = copy_to(obj, dir_path, name)
    except FileNotFoundError as e:
        new_obj = None

    return new_obj
