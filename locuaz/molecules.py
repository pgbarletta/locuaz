from itertools import chain
from functools import singledispatch
import logging
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

    file: FileHandle = field()
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

    # TODO: DEPRECATE
    # def update_chains(self, *, target_chains: List, binder_chains: List):
    #     self.target_chains = {chainID: None for chainID in target_chains}
    #     self.binder_chains = {chainID: None for chainID in binder_chains}

    # @classmethod
    # def from_config(cls, config: Dict):
    #     top_path = Path(config["paths"]["input"], config["md"]["topology"])
    #     self = cls.from_path(top_path)
    #     return self

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
    # TODO: DEPRECATE
    # @classmethod
    # def from_config(cls, config: Dict) -> "ZipTopology":
    #     self = super().from_config(config)
    #     self.target_chains = tuple(config["target"]["chainID"])
    #     self.binder_chains = tuple(config["binder"]["chainID"])

    #     return self


@define
class TrrTrajectory(Trajectory):
    pass


@define
class XtcTrajectory(Trajectory):
    pass


# TODO: deprecate
@define
class GROTopology(Topology):
    dir_path: Path = field(init=False)
    # TODO: This overrides topology's attributes?
    target_chains_files: Dict[str, FileHandle] = field(init=False)
    binder_chains_files: Dict[str, FileHandle] = field(init=False)

    def add_chain_tops(
        self,
        target_chainid: List[str],
        target_top: List[str],
        binder_chainid: List[str],
        binder_top: List[str],
    ):
        self.dir_path = self.file.path.parent
        self.target_chains_files = dict()
        for chainID, top in zip(target_chainid, target_top):
            self.target_chains_files[chainID] = FileHandle(self.dir_path / top)

        self.binder_chains_files = dict()
        for chainID, top in zip(binder_chainid, binder_top):
            self.binder_chains_files[chainID] = FileHandle(self.dir_path / top)

    # TODO: DEPRECATE
    # @classmethod
    # def from_config(cls, config: Dict) -> "GROTopology":
    #     # top_path = Path(config["paths"]["input"], config["md"]["topology"])
    #     # self = cls.from_path(top_path)
    #     self = super().from_config(config)
    #     self.add_chain_tops(
    #         config["target"]["chainID"],
    #         config["md"]["target_topology"],
    #         config["binder"]["chainID"],
    #         config["md"]["binder_topology"],
    #     )
    #     return self

    @classmethod
    def from_topology(
        cls, top_file: FileHandle, target_chains_files: Dict, binder_chains_files: Dict
    ):
        self = cls(top_file)
        self.target_chains_files = target_chains_files
        self.binder_chains_files = binder_chains_files

        return self


@define
class AbstractComplex(metaclass=ABCMeta):
    name: str = field(converter=str)
    dir_handle: DirHandle = field(converter=DirHandle)  # type: ignore
    pdb: PDBStructure = field()
    top: Topology = field()
    tra: Trajectory = field(init=False, default=None)

    
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

    # TODO: `from_gro_zip`'s is GROMACS specific.
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
        gmx_bin: str = "gmx",
    ) -> "AbstractComplex":
        raise NotImplementedError

    def __str__(self) -> str:
        return str(self.dir_handle)


# TODO: make frozen, put all attributes as init = true and kw_only
@define
class GROComplex(AbstractComplex):
    gro: GROStructure = field()
    tpr: TPRFile = field(init=False, default=None)
    ndx: FileHandle = field(init=False, default=None)

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
            logging.info(
                f"Creating GROComplex from: {in_pdb}. This PDB will be"
                "backed up and replaced with a fixed up PDB."
            )

            temp_pdb = PDBStructure.from_path(in_pdb)
        except Exception as e:
            print(
                f"Could not get input PDB file from: {str(input_dir)}",
                flush=True,
            )
            raise e
        else:
            try:
                str_pdb, str_gro, top = get_gro_ziptop_from_pdb(
                    pdb=temp_pdb,
                    target_chains=target_chains,
                    binder_chains=binder_chains,
                    gmx_path=gmx_bin,
                    add_ions=True,
                    water_type=water_type,
                    force_field=force_field,
                )
                tpr = get_tpr(gro=str_gro, top=top, gmx_bin=gmx_bin)

            except Exception as e:
                print(
                    f"Could not get zip tology .gro or TPR files from: "
                    f"{temp_pdb.file.path}",
                    flush=True,
                )
                raise e
            else:
                self = GROComplex(name, input_dir, str_pdb, top, str_gro)
                self.tpr = tpr
                self.update_all_ndxs(gmx_bin=gmx_bin)

                return self

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
                f"Could not get input .gro and .zip files from: {str(input_dir)}",
                flush=True,
            )
            raise e
        else:
            try:
                pdb, tpr = get_pdb_tpr(gro=gro, top=top, gmx_bin=gmx_bin)

            except Exception as e:
                print(
                    f"Could not get PDB and TPR files from: "
                    f"{gro.file.path} and {top.file.path}",
                    flush=True,
                )
                raise e
            else:
                self = GROComplex(name, input_dir, pdb, top, gro)
                self.tpr = tpr
                self.update_all_ndxs(gmx_bin=gmx_bin)

                return self

    @classmethod
    def from_complex(
        cls,
        *,
        name: str,
        iter_path: Path,
        target_chains: Sequence,
        binder_chains: Sequence,
        gmx_bin: str = "gmx",
    ) -> "GROComplex":
        try:
            str_pdb = PDBStructure.from_path(iter_path / (name + ".pdb"))
            str_gro = GROStructure.from_path(iter_path / (name + ".gro"))
            top = ZipTopology.from_path(iter_path / (name + ".zip"))
            top.target_chains = tuple(target_chains)
            top.binder_chains = tuple(binder_chains)
            try:
                traj = XtcTrajectory.from_path(iter_path / (name + ".xtc"))
            except FileNotFoundError as e:
                traj = TrrTrajectory.from_path(iter_path / (name + ".trr"))
            tpr = TPRFile.from_path(iter_path / (name + ".tpr"))
        except Exception as e:
            logging.error(
                f"Could not get some file from: {iter_path}\n"
                f"Input parameters:"
                f"\n\t{name}\n\t{iter_path}"
                f"\n\t{target_chains}\n\t{binder_chains}"
            )
            raise e
        else:
            self = GROComplex(name, iter_path, str_pdb, top, str_gro)
            self.tra = traj
            self.tpr = tpr
            self.update_all_ndxs(gmx_bin=gmx_bin)

            return self

    def update_all_ndxs(self, gmx_bin: str = "gmx"):
        try:
            # First, get the target and binder chains
            target_ndx = self.write_target_ndx()
            binder_ndx = self.write_binder_ndx()

            # Then, cat them into 1 file
            target_and_binder = catenate(
                self.dir_handle.dir_path / "target_and_binder.ndx",
                target_ndx,
                binder_ndx,
            )

            # Use `target_and_binder` to get a selection of both
            target_and_binder_both = (
                Path(self.dir_handle) / "target_and_binder_both.ndx"
            )
            selector = MakeNdx(
                input_structure_path=str(self.pdb.file.path),
                input_ndx_path=str(target_and_binder.path),
                output_ndx_path=str(target_and_binder_both),
                properties={"selection": '"target" | "binder"'},
            )
            launch_biobb(selector)

            # Finally, also add the negation of target+binder
            complex_ndx_fn = Path(self.dir_handle) / f"{self.name}.ndx"
            selector = MakeNdx(
                input_structure_path=str(self.pdb.file.path),
                input_ndx_path=str(target_and_binder_both),
                output_ndx_path=str(complex_ndx_fn),
                properties={"selection": '! "target_binder"'},
            )
            launch_biobb(selector)

            # Save the final ndx file.
            self.ndx = FileHandle(complex_ndx_fn)

            # Clean-up
            target_ndx.unlink()
            binder_ndx.unlink()
            target_and_binder.unlink()
            target_and_binder_both.unlink()

        except Exception as e:
            raise e

    def write_target_ndx(self, *, selname: str = "target") -> FileHandle:
        selection = " or ".join(
            [f"chain {chainID}" for chainID in self.top.target_chains]
        )
        ndx_fn = Path(self.dir_handle) / f"{selname}.ndx"

        selector = Gmxselect(
            input_structure_path=str(self.pdb.file.path),
            output_ndx_path=str(ndx_fn),
            properties={"selection": selection},
        )
        launch_biobb(selector)

        ndx = FileHandle(ndx_fn)
        update_header(ndx, f"[ {selname} ]\n")

        return ndx

    def write_binder_ndx(self, *, selname: str = "binder") -> FileHandle:
        selection = " or ".join(
            [f"chain {chainID}" for chainID in self.top.binder_chains]
        )
        ndx_fn = Path(self.dir_handle) / f"{selname}.ndx"

        selector = Gmxselect(
            input_structure_path=str(self.pdb.file.path),
            output_ndx_path=str(ndx_fn),
            properties={"selection": selection},
        )
        launch_biobb(selector)

        ndx = FileHandle(ndx_fn)
        update_header(ndx, f"[ {selname} ]\n")

        return ndx


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
    """get_gro_ziptop_from_pdb() does a pdb2gmx from the PDB and tries to keep
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
        "dev": "-renum",
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

    # Get PDB from GRO file. Gromacs should back up the older input PDB
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
    nonwat_pdb_fn = Path(complex.dir_handle) / ("nonwat_" + complex.name + ".pdb")
    get_protein = GMXTrjConvStr(
        input_structure_path=str(complex.pdb.file.path),
        input_top_path=str(complex.tpr.file.path),
        input_index_path=str(complex.ndx.path),
        output_str_path=str(nonwat_pdb_fn),
        properties={"selection": "target_binder"},
    )
    launch_biobb(get_protein)
    nonwat_pdb = PDBStructure.from_path(nonwat_pdb_fn)

    # Water and ions
    wation_pdb_fn = Path(complex.dir_handle) / ("wation_" + complex.name + ".pdb")
    get_water_ions = GMXTrjConvStr(
        input_structure_path=str(complex.pdb.file.path),
        input_top_path=str(complex.tpr.file.path),
        input_index_path=str(complex.ndx.path),
        output_str_path=str(wation_pdb_fn),
        properties={"selection": "!target_binder"},
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


def catenate_pdbs(*pdbs: PDBStructure, pdb_out_path: Path, gmx_bin: str = "gmx"):
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


@singledispatch
def copy_mol_to(obj, dir_path: Path, name=None):
    raise NotImplementedError


@copy_mol_to.register
def _(obj: Structure, dir_path: Path, name=None):
    new_file = copy_to(obj.file, dir_path, name=name)
    return Structure(new_file)


@copy_mol_to.register
def _(obj: GROTopology, dir_path: Path, name=None):
    new_file = copy_to(obj.file, dir_path, name=name)
    target_chains = {
        chainID: copy_to(chain_file, dir_path)
        for chainID, chain_file in obj.target_chains_files.items()
    }
    binder_chains = {
        chainID: copy_to(chain_file, dir_path)
        for chainID, chain_file in obj.binder_chains_files.items()
    }

    return GROTopology.from_topology(new_file, target_chains, binder_chains)


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
    str_gro = copy_to(obj.gro, dir_path, name)
    top = copy_to(obj.top, dir_path, name)

    new_cpx = GROComplex(obj.name, dir_path, str_pdb, top, str_gro)

    if obj.tra:
        traj = copy_to(obj.tra, dir_path, name)
        GROComplex.tra = traj
    if obj.ndx:
        ndx = copy_to(obj.ndx, dir_path, name)
        GROComplex.ndx = ndx

    return new_cpx
