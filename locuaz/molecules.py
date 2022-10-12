from itertools import chain
from functools import singledispatch
import MDAnalysis as mda
from pathlib import Path

from attrs import define, field
from fileutils import DirHandle, FileHandle, copy_to, catenate, update_header
from abc import ABCMeta, abstractmethod
from typing import List, Sequence, Set, Dict, Tuple, Optional
from biobb_md.gromacs.pdb2gmx import Pdb2gmx
from biobb_md.gromacs.gmxselect import Gmxselect
from biobb_md.gromacs.make_ndx import MakeNdx
from biobb_md.gromacs.grompp import Grompp
from biobb_analysis.gromacs.gmx_trjconv_str import GMXTrjConvStr
from biobb_md.gromacs.genion import Genion
from primitives import launch_biobb


# fmt: off
res_1 = ("G", "A", "L", "M", "F", "W", "K", "Q", "E", "S",
    "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T")
d_res_1 = dict(zip(res_1, range(0, len(res_1))))

res_3 = ("Gly", "Ala", "Leu", "Met", "Phe", "Trp", "Lys", "Gln", "Glu",
    "Ser", "Pro", "Val", "Ile", "Cys", "Tyr", "His", "Arg", "Asn", "Asp", "Thr")
d_res_3 = dict(zip(res_3, range(0, len(res_3))))
# fmt: on


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
    def from_path(cls, path: Path):
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

    @classmethod
    def from_path_with_chains(
        cls, path: Path, *, target_chains: Sequence, binder_chains: Sequence
    ) -> "Topology":
        self = super().from_path(path)
        self.target_chains = tuple(target_chains)
        self.binder_chains = tuple(binder_chains)
        return self


@define
class AmberTopology(Topology):
    pass


@define
class ZipTopology(Topology):
    pass


@define
class TrrTrajectory(Trajectory):
    pass


@define
class XtcTrajectory(Trajectory):
    pass


@define
class AbstractComplex(metaclass=ABCMeta):
    name: str = field(converter=str)
    dir: DirHandle = field(converter=DirHandle, kw_only=True)  # type: ignore
    pdb: PDBStructure = field(kw_only=True)
    top: Topology = field(kw_only=True, default=None)
    tra: Trajectory = field(kw_only=False, default=None)

    @classmethod
    @abstractmethod
    def from_pdb(
        cls,
        *,
        name: str,
        input_dir: Path,
        target_chains: Sequence,
        binder_chains: Sequence,
        gmx_bin: str = "gmx",
        water_type: str = "tip3p",
        force_field: str = "amber99sb-ildn",
    ) -> "AbstractComplex":
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def from_gro_zip(
        cls,
        *,
        name: str,
        input_dir: Path,
        target_chains: Sequence,
        binder_chains: Sequence,
        gmx_bin: str = "gmx",
    ) -> "AbstractComplex":
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def from_complex(
        cls,
        *,
        name: str,
        iter_path: Path,
        target_chains: Sequence,
        binder_chains: Sequence,
        ignore_cpt: bool = True,
        gmx_bin: str = "gmx",
    ) -> "AbstractComplex":
        raise NotImplementedError

    def get_cryst1_record(self) -> str:
        return self.pdb.get_cryst1_record()

    def __str__(self) -> str:
        return str(self.dir)


@define(frozen=True)
class GROComplex(AbstractComplex):
    gro: GROStructure = field(kw_only=True, default=None)
    tpr: TPRFile = field(kw_only=False, default=None)
    cpt: FileHandle = field(kw_only=False, default=None)
    ndx: FileHandle = field(kw_only=False, default=None)

    @classmethod
    def from_pdb(
        cls,
        *,
        name: str,
        input_dir: Path,
        target_chains: Sequence,
        binder_chains: Sequence,
        gmx_bin: str = "gmx",
        water_type: str = "tip3p",
        force_field: str = "amber99sb-ildn",
    ) -> "GROComplex":
        try:
            in_pdb = input_dir / (name + ".pdb")
            # This PDB will be backed up and replaced with a fixed up PDB.
            temp_pdb = PDBStructure.from_path(in_pdb)
        except Exception as e:
            print(f"Could not get input PDB file from: {input_dir}", flush=True)
            raise e
        try:
            # The PDB should have box info in it, so it'll be pasted onto the GRO. Else,
            # pdb2gmx will compute the binding box of the PDB and use that one.
            pdb, gro, top = get_gro_ziptop_from_pdb(
                pdb=temp_pdb,
                target_chains=target_chains,
                binder_chains=binder_chains,
                gmx_path=gmx_bin,
                add_ions=True,
                water_type=water_type,
                force_field=force_field,
            )
        except Exception as e:
            print(f"Could not get zip tology .gro files from: {temp_pdb}", flush=True)
            raise e
        try:
            traj = XtcTrajectory.from_path(input_dir / (name + ".xtc"))
        except FileNotFoundError as e:
            try:
                traj = TrrTrajectory.from_path(input_dir / (name + ".trr"))
            except FileNotFoundError as ee:
                traj = None
        tpr = get_tpr(gro=gro, top=top, gmx_bin=gmx_bin)
        ndx = generate_ndx(
            name,
            pdb=pdb,
            target_chains=target_chains,
            binder_chains=binder_chains,
        )

        return GROComplex(
            name,
            dir=input_dir,
            pdb=pdb,
            top=top,
            tra=traj,
            gro=gro,
            tpr=tpr,
            ndx=ndx,
        )

    @classmethod
    def from_gro_zip(
        cls,
        *,
        name: str,
        input_dir: Path,
        target_chains: Sequence,
        binder_chains: Sequence,
        gmx_bin: str = "gmx",
    ) -> "AbstractComplex":
        try:
            gro = GROStructure.from_path(input_dir / (name + ".gro"))

            top = ZipTopology.from_path(input_dir / (name + ".zip"))
            top.target_chains = tuple(target_chains)
            top.binder_chains = tuple(binder_chains)
        except Exception as e:
            print(
                f"Could not get input .gro and .zip files from: {input_dir}", flush=True
            )
            raise e

        pdb, tpr = get_pdb_tpr(gro=gro, top=top, gmx_bin=gmx_bin)
        try:
            traj = XtcTrajectory.from_path(input_dir / (name + ".xtc"))
        except FileNotFoundError as e:
            try:
                traj = TrrTrajectory.from_path(input_dir / (name + ".trr"))
            except FileNotFoundError as ee:
                traj = None
        ndx = generate_ndx(
            name,
            pdb=pdb,
            target_chains=target_chains,
            binder_chains=binder_chains,
        )

        return GROComplex(
            name,
            dir=input_dir,
            pdb=pdb,
            top=top,
            tra=traj,
            gro=gro,
            tpr=tpr,
            ndx=ndx,
        )

    @classmethod
    def from_complex(
        cls,
        *,
        name: str,
        iter_path: Path,
        target_chains: Sequence,
        binder_chains: Sequence,
        ignore_cpt: bool = True,
        gmx_bin: str = "gmx",
    ) -> "GROComplex":
        try:
            pdb = PDBStructure.from_path(iter_path / (name + ".pdb"))
            gro = GROStructure.from_path(iter_path / (name + ".gro"))
            top = ZipTopology.from_path(iter_path / (name + ".zip"))
            top.target_chains = tuple(target_chains)
            top.binder_chains = tuple(binder_chains)
            try:
                traj = XtcTrajectory.from_path(iter_path / (name + ".xtc"))
            except FileNotFoundError:
                try:
                    traj = TrrTrajectory.from_path(iter_path / (name + ".trr"))
                except FileNotFoundError as ee:
                    traj = None
            tpr = TPRFile.from_path(iter_path / (name + ".tpr"))
            if ignore_cpt:
                cpt = None
            else:
                try:
                    cpt = FileHandle(iter_path / (name + ".cpt"))
                except FileNotFoundError as e:
                    cpt = None
        except Exception as ee:
            raise RuntimeError("from_complex() failed.") from ee
        ndx = generate_ndx(
            name,
            pdb=pdb,
            target_chains=target_chains,
            binder_chains=binder_chains,
        )

        return GROComplex(
            name,
            dir=iter_path,
            pdb=pdb,
            top=top,
            tra=traj,
            gro=gro,
            tpr=tpr,
            cpt=cpt,
            ndx=ndx,
        )


def get_gro_ziptop_from_pdb(
    *,
    pdb: PDBStructure,
    target_chains: Sequence,
    binder_chains: Sequence,
    gmx_path: str = "gmx",
    water_type: str = "tip3p",
    force_field: str = "amber99sb-ildn",
    add_ions: bool = False,
) -> Tuple[PDBStructure, GROStructure, ZipTopology]:
    """get_gro_ziptop_from_pdb does a pdb2gmx from the PDB and tries to keep
    the system neutral, which may alter the topology so a new PDB will be
    written with the same name as the original, which will be backed up by GROMACS.

    Args:
        pdb (PDBStructure): input PDB
        target_chains (Sequence): these will be used to construct the ZipTopology
        binder_chains (Sequence): these will be used to construct the ZipTopology
        gmx_path (str, optional): for all the biobb tools. Defaults to "gmx".
        water_type (str, optional): argument to pdb2gmx. Defaults to "tip3p".
        force_field (str, optional): argument to pdb2gmx. Defaults to "amber99sb-ildn".

    Returns:
        Tuple[PDBStructure, GROStructure, ZipTopology]: Proper, nice, system.
    """

    local_dir = pdb.file.path.parent
    name = pdb.name

    # Generate the first set of GROStructure and ZipTopology
    props = {
        "gmx_path": gmx_path,
        "water_type": water_type,
        "force_field": force_field,
        "ignh": True,
        # "dev": "-renum",
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

    if add_ions:
        # Build a .tpr for genion
        gen_tpr_fn = local_dir / ("genion_" + name + ".tpr")
        grompepe = Grompp(
            input_gro_path=str(pre_gro_fn),
            input_top_zip_path=str(pre_top_fn),
            output_tpr_path=str(gen_tpr_fn),
            properties={"gmx_path": gmx_path, "maxwarn": 2},
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
            properties={"gmx_path": gmx_path, "neutral": True, "concentration": 0.0},
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
            "gmx_path": gmx_path,
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
        properties={"gmx_path": gmx_path},
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
        properties={"gmx_path": gmx_bin, "simulation_type": "index"},
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
        properties={"gmx_path": gmx_bin},
    )
    launch_biobb(trjconv)

    return PDBStructure.from_path(pdb_fn), tpr


def split_solute_and_solvent(complex: GROComplex) -> Tuple[PDBStructure, PDBStructure]:
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
        properties={"selection": "complex"},
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
        properties={"selection": "solvent_ions"},
    )
    launch_biobb(get_water_ions)

    wation_pdb = PDBStructure.from_path(wation_pdb_fn)

    return nonwat_pdb, wation_pdb


def fix_pdb(
    pdb_in: PDBStructure, pdb_out_path: Path, gmx_bin: str = "gmx"
) -> PDBStructure:
    # This is quite cumbersome.
    wrk_dir = pdb_in.file.path.parent

    # First, use pdb2gmx to fix the PDB file. There's no way of ignoring
    # Hydrogens and getting a good output PDB, so the output has to be
    # a .gro file
    props = {
        "gmx_path": str(gmx_bin),
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
        properties={"gmx_path": gmx_bin, "maxwarn": 1},
    )
    launch_biobb(grompepe)

    # Finally, get a PDB file from the GRO file
    trjconv = GMXTrjConvStr(
        input_structure_path=str(temp_gro),
        input_top_path=str(temp_tpr),
        output_str_path=str(pdb_out_path),
        properties={"gmx_path": gmx_bin},
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
    *, pdb: PDBStructure, chains: Sequence, selname: str, gmx_bin: str = "gmx"
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
        uni_pdb = mda.Universe(pdb)
        ndx_file = Path(Path(pdb).parent, f"{name}.ndx")

        uni_pdb.select_atoms(
            " or ".join([f"segid {chainID}" for chainID in target_chains])
        ).write(ndx_file, name="target", mode="w")
        uni_pdb.select_atoms(
            " or ".join([f"segid {chainID}" for chainID in binder_chains])
        ).write(ndx_file, name="binder", mode="a")
        # TODO: this won't work with glyco mods and stuff
        uni_pdb.select_atoms("protein").write(ndx_file, name="complex", mode="a")
        uni_pdb.select_atoms("not protein").write(
            ndx_file, name="solvent_ions", mode="a"
        )
        uni_pdb.select_atoms("all").write(ndx_file, name="sistema", mode="a")

        return FileHandle(ndx_file)

    except Exception as e:
        raise e


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


@copy_mol_to.register
def _(obj: GROComplex, dir_path: Path, name=None):
    str_pdb = copy_to(obj.pdb, dir_path, name)
    top = try_copy_to(obj.top, dir_path, name)
    tra = try_copy_to(obj.tra, dir_path, name)
    gro = try_copy_to(obj.gro, dir_path, name)
    tpr = try_copy_to(obj.tpr, dir_path, name)
    ndx = try_copy_to(obj.ndx, dir_path, name)

    return GROComplex(
        obj.name, dir=dir_path, pdb=str_pdb, top=top, tra=tra, gro=gro, tpr=tpr, ndx=ndx
    )


def try_copy_to(obj, dir_path: Path, name=None):
    try:
        new_obj = copy_to(obj, dir_path, name)
    except FileNotFoundError as e:
        new_obj = None

    return new_obj
