from pathlib import Path
from abc import ABCMeta
from abc import ABCMeta, abstractmethod
from typing import Sequence, Dict, Optional, Union

from attrs import define, field

from fileutils import DirHandle, FileHandle, copy_to
from molecules import (
    PDBStructure,
    GROStructure,
    Topology,
    ZipTopology,
    TPRFile,
    Trajectory,
    XtcTrajectory,
    TrrTrajectory,
    get_gro_ziptop_from_pdb,
    generate_ndx,
    get_tpr,
    get_pdb_tpr,
    copy_mol_to,
    try_copy_to,
)


@define(frozen=True)
class AbstractComplex(metaclass=ABCMeta):
    name: str = field(converter=str)
    dir: Union[DirHandle, Path] = field(converter=DirHandle, kw_only=True)  # type: ignore
    pdb: PDBStructure = field(kw_only=True)
    top: Topology = field(kw_only=True, default=None)
    tra: Optional[Trajectory] = field(kw_only=False, default=None)

    @classmethod
    @abstractmethod
    def from_pdb(
        cls,
        *,
        name: str,
        input_dir: Path,
        target_chains: Sequence,
        binder_chains: Sequence,
        md_config: Dict,
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
    """GROComplex Main model of the optimized complex. It's immutable and always ready for
    MD. Each time a Complex is built (either from a PDB, a set of GRO and ZIP topology,
    or another complex), all its attributes are filled out and its associated files are
    created. Some of them, as checkpoint and trajectoy, are optional.

    Attributes:
        name (str): name of the complex. All of the associated files (through filehandles)
            will be named as f'name.{extension}'.
        dir (Union[DirHandle, Path]): path to the complex.
        pdb (PDBStructure):
        top (ZipTopology): zip topolo
        tra (Optional[Trajectory]):
        gro (GROStructure):
        tpr (TPRFile):
        cpt (Optional[FileHandle]):
        ndx (FileHandle): index file with selections for target, binder, protein,
            not protein and System.
    Raises:
        RuntimeError: failure to build the complex.
    """

    gro: GROStructure = field(kw_only=True, default=None)
    tpr: TPRFile = field(kw_only=False, default=None)
    cpt: Optional[FileHandle] = field(kw_only=False, default=None)
    ndx: FileHandle = field(kw_only=False, default=None)

    @classmethod
    def from_pdb(
        cls,
        *,
        name: str,
        input_dir: Union[Path, DirHandle],
        target_chains: Sequence,
        binder_chains: Sequence,
        md_config: Dict,
    ) -> "GROComplex":
        gmx_bin: str = md_config.get("gmx_bin", "gmx")
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
                md_config=md_config,
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
    ) -> "GROComplex":
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
