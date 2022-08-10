from pathlib import Path
from attrs import define, field
from fileutils import DirHandle, FileHandle
from abc import ABCMeta, abstractmethod
from typing import List, Set, Dict, Tuple, Optional
from biobb_md.gromacs.pdb2gmx import pdb2gmx

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


class Topology(AbstractFileObject):
    target_chains: Dict[str, None] = field(init=False)
    binder_chains: Dict[str, None] = field(init=False)

    def update_chains(self, *, target_chains: List, binder_chains: List):
        self.target_chains = {chainID: None for chainID in target_chains}
        self.binder_chains = {chainID: None for chainID in binder_chains}

    @classmethod
    def from_config(cls, config: Dict):
        top_path = Path(config["paths"]["input"], config["md"]["topology"])
        self = cls.from_path(top_path)
        return self

    @classmethod
    def from_path_with_chains(
        cls, path: Path, *, target_chains: List, binder_chains: List
    ):
        # file = FileHandle(path)
        # self = cls(file)
        self = super().from_path(path)
        self.update_chains(target_chains=target_chains, binder_chains=binder_chains)
        return self


@define
class AmberTopology(Topology):
    pass


@define
class ZipTopology(Topology):
    @classmethod
    def from_config(cls, config: Dict) -> "ZipTopology":
        self = super().from_config(config)
        self.update_chains(
            target_chains=config["target"]["chainID"],
            binder_chains=config["binder"]["chainID"],
        )
        return self


@define
class TrrTrajectory(Trajectory):
    pass


@define
class XtcTrajectory(Trajectory):
    pass


@define
class GROTopology(Topology):
    dir_path: Path = field(init=False)
    # TODO: This overrides topology's attributes?
    target_chains: Dict[str, FileHandle] = field(init=False)
    binder_chains: Dict[str, FileHandle] = field(init=False)

    def add_chain_tops(
        self,
        target_chainid: List[str],
        target_top: List[str],
        binder_chainid: List[str],
        binder_top: List[str],
    ):
        self.dir_path = self.file.path.parent
        self.target_chains = dict()
        for chainID, top in zip(target_chainid, target_top):
            self.target_chains[chainID] = FileHandle(self.dir_path / top)

        self.binder_chains = dict()
        for chainID, top in zip(binder_chainid, binder_top):
            self.binder_chains[chainID] = FileHandle(self.dir_path / top)

    @classmethod
    def from_config(cls, config: Dict) -> "GROTopology":
        # top_path = Path(config["paths"]["input"], config["md"]["topology"])
        # self = cls.from_path(top_path)
        self = super().from_config(config)
        self.add_chain_tops(
            config["target"]["chainID"],
            config["md"]["target_topology"],
            config["binder"]["chainID"],
            config["md"]["binder_topology"],
        )
        return self

    @classmethod
    def from_topology(
        cls, top_file: FileHandle, target_chains: Dict, binder_chains: Dict
    ):
        self = cls(top_file)
        self.target_chains = target_chains
        self.binder_chains = binder_chains

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
    def from_config(cls, config: Dict) -> "AbstractComplex":
        raise NotImplementedError


@define
class GROComplex(AbstractComplex):
    gro: GROStructure = field()
    tpr: TPRFile = field(init=False, default=None)
    target_ndx: FileHandle = field(init=False, default=None)
    binder_ndx: FileHandle = field(init=False, default=None)
    complex_ndx: FileHandle = field(init=False, default=None)

    @classmethod
    def from_config(cls, config: Dict) -> "GROComplex":
        try:
            dir_path = Path(config["paths"]["input"])
            str_pdb = PDBStructure.from_path(
                dir_path / (config["main"]["name"] + ".pdb")
            )
            str_gro = GROStructure.from_path(
                dir_path / (config["main"]["name"] + ".gro")
            )

            top: Topology
            if Path(config["md"]["topology"]).suffix == ".zip":
                top = ZipTopology.from_config(config)
            else:
                top = GROTopology.from_config(config)
        except Exception as e:
            print(
                f'Could not get input files from: {config["paths"]["input"]}',
                flush=True,
            )
            raise e
        else:
            return GROComplex(config["main"]["name"], dir_path, str_pdb, top, str_gro)

    @classmethod
    def from_iter(
        cls, name: str, iter_path: Path, target_chains: List, binder_chains: List
    ) -> "GROComplex":
        try:
            str_pdb = PDBStructure.from_path(iter_path / ("npt_" + name + ".pdb"))
            str_gro = GROStructure.from_path(iter_path / ("npt_" + name + ".gro"))
            top = ZipTopology.from_path(iter_path / ("wet_topol.zip"))
            top.update_chains(target_chains=target_chains, binder_chains=binder_chains)
            traj = XtcTrajectory.from_path(iter_path / ("npt_" + name + ".xtc"))
            tpr = TPRFile.from_path(iter_path / ("npt_" + name + ".tpr"))
        except Exception as e:
            print(
                f"Could not get input files from: {iter_path}",
                flush=True,
            )
            raise e
        else:
            self = GROComplex(name, iter_path, str_pdb, top, str_gro)
            self.tra = traj
            self.tpr = tpr
            return self


def get_zip_topology(
    *,
    pdb: PDBStructure,
    gmx_path: str = "gmx",
    water_type: str = "tip3p",
    force_field: str = "amber99sb-ildn",
) -> Tuple[GROStructure, ZipTopology]:

    props = {
        "gmx_path": gmx_path,
        "water_type": water_type,
        "force_field": force_field,
        "ignh": True,
    }
    gro = pdb.file.path.parent / (pdb.name + ".gro")
    top = pdb.file.path.parent / (pdb.name + ".zip")
    pdb2gmx(
        input_pdb_path=str(pdb.file),
        output_gro_path=str(gro),
        output_top_zip_path=str(top),
        properties=props,
    )

    return GROStructure.from_path(gro), ZipTopology.from_path(top)
