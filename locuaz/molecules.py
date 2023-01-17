from itertools import chain
from functools import singledispatch
import MDAnalysis as mda
from pathlib import Path
import numpy as np
from abc import ABCMeta
from collections.abc import Iterable
from typing import List, Dict, Tuple, Any

from attrs import define, field

from biobb_gromacs.gromacs.pdb2gmx import Pdb2gmx
from biobb_gromacs.gromacs.gmxselect import Gmxselect
from biobb_gromacs.gromacs.grompp import Grompp
from biobb_analysis.gromacs.gmx_trjconv_str import GMXTrjConvStr


from fileutils import FileHandle, copy_to, update_header
from primitives import launch_biobb


@define
class AbstractFileObject(metaclass=ABCMeta):
    """AbstractFileObject Abstract Base Class for all objects associated to a file

    Args:
        path ([str, Path]): path to file, will use to construct the FileHandle object
    Attributes:
        file (FileHandle): FileHandle to the object
        name (str): name of the file
        extension (str): extension, filetype, of the file
    """

    file: FileHandle = field(converter=FileHandle)  # type: ignore
    name: str = field(init=False)
    ext: str = field(init=False)

    # def __init__(self, file_str: Path):
    #     file_path = Path(file_str)
    #     self.file = FileHandle(file_path)
    #     self.name, self.ext = file_path.name.split(".")

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
        name (str):
        file (FileHandle): path to the file
    """

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
        name (str):
        file (FileHandle): path to the file
    """

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
    selection_protein: str = field(init=False)
    selection_not_protein: str = field(init=False)

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
        self.selection_protein = " or ".join(
            ["protein", self.selection_target, self.selection_binder]
        )
        self.selection_not_protein = " and not ".join(
            ["not protein", self.selection_target, self.selection_binder]
        )

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
        self.selection_protein = " or ".join(
            ["protein", self.selection_target, self.selection_binder]
        )
        self.selection_not_protein = " and not ".join(
            ["not protein", self.selection_target, self.selection_binder]
        )

        return self


@define
class TrrTrajectory(Trajectory):
    pass


@define
class XtcTrajectory(Trajectory):
    pass


def get_gro_ziptop_from_pdb(
    *,
    pdb: PDBStructure,
    target_chains: Sequence,
    binder_chains: Sequence,
    md_config: Dict,
    add_ions: bool = False,
) -> Tuple[PDBStructure, GROStructure, ZipTopology]:
    """get_gro_ziptop_from_pdb does a pdb2gmx from the PDB and tries to keep
    the system neutral, which may alter the topology so a new PDB will be
    written with the same name as the original, which will be backed up by GROMACS.

    Args:
        pdb (PDBStructure): input PDB
        target_chains (Sequence): these will be used to construct the ZipTopology
        binder_chains (Sequence): these will be used to construct the ZipTopology
        binary_path (str, optional): for all the biobb tools. Defaults to "gmx".
        water_type (str, optional): argument to pdb2gmx. Defaults to "tip3p".
        force_field (str, optional): argument to pdb2gmx. Defaults to "amber99sb-ildn".

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
        gen_tpr_fn.unlink()
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

    top = ZipTopology.from_path(top_fn)
    top.target_chains = tuple(target_chains)
    top.binder_chains = tuple(binder_chains)

    return PDBStructure.from_path(pdb_fn), GROStructure.from_path(gro_fn), top


def get_tpr(
    *,
    gro: GROStructure,
    top: ZipTopology,
    gmx_bin: str = "gmx",
) -> TPRFile:

    tpr = Path(gro.file.path.parent) / (gro.name + ".tpr")
    grompepe = Grompp(
        input_gro_path=str(gro.file.path),
        input_top_zip_path=str(top.file.path),
        output_tpr_path=str(tpr),
        properties={"binary_path": gmx_bin, "simulation_type": "index"},
    )
    launch_biobb(grompepe)

    return TPRFile.from_path(tpr)


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


def fix_pdb(
    pdb_in: PDBStructure, pdb_out_path: Path, gmx_bin: str = "gmx"
) -> PDBStructure:
    # This is quite cumbersome.
    wrk_dir = pdb_in.file.path.parent

    # First, use pdb2gmx to fix the PDB file. There's no way of ignoring
    # Hydrogens and getting a good output PDB, so the output has to be
    # a .gro file
    props = {
        "binary_path": str(gmx_bin),
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
    out_pdb = fix_pdb(PDBStructure.from_path(temp_pdb), pdb_out_path, gmx_bin)
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
    target_chains: Sequence,
    binder_chains: Sequence,
) -> FileHandle:
    try:
        uni_pdb = mda.Universe(str(pdb))
        ndx_file = Path(Path(pdb).parent, f"{name}.ndx")

        uni_pdb.select_atoms(
            " or ".join([f"segid {chainID}" for chainID in target_chains])
        ).write(ndx_file, name="target", mode="w")
        uni_pdb.select_atoms(
            " or ".join([f"segid {chainID}" for chainID in binder_chains])
        ).write(ndx_file, name="binder", mode="a")
        # TODO: this won't work with glyco mods and stuff
        uni_pdb.select_atoms("protein").write(ndx_file, name="Protein", mode="a")
        uni_pdb.select_atoms("not protein").write(
            ndx_file, name="Non-Protein", mode="a"
        )
        uni_pdb.select_atoms("all").write(ndx_file, name="sistema", mode="a")

        return FileHandle(ndx_file)

    except Exception as e:
        raise e


def read_ndx(ndx_path: Path) -> Dict:
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
def copy_mol_to(obj, dir_path: Path, name=None):
    raise NotImplementedError


@copy_mol_to.register
def _(obj: Structure, dir_path: Path, name=None):
    new_file = copy_to(obj.file, dir_path, name=name)
    return Structure(new_file)


@copy_mol_to.register
def _(obj: ZipTopology, dir_path: Path, name=None):
    new_file = copy_to(obj.file, dir_path, name=name)
    ZipTopology.from_path_with_chains(
        new_file.path, target_chains=obj.target_chains, binder_chains=obj.binder_chains
    )


@copy_mol_to.register
def _(obj: Trajectory, dir_path: Path, name=None):
    new_file = copy_to(obj.file, dir_path, name=name)
    return Trajectory(new_file)


def try_copy_to(obj, dir_path: Path, name=None) -> Any:
    try:
        new_obj = copy_to(obj, dir_path, name)
    except FileNotFoundError as e:
        new_obj = None

    return new_obj
