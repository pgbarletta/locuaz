from abc import ABCMeta, abstractmethod
from collections.abc import Iterable
from pathlib import Path
from typing import Dict, Optional, Union

from attrs import define, field

from locuaz.fileutils import DirHandle, FileHandle, copy_to
from locuaz.molecules import (
    PDBStructure,
    GROStructure,
    Topology,
    ZipTopology,
    TPRFile,
    Trajectory,
    XtcTrajectory,
    TrrTrajectory,
    generate_ndx,
    get_tpr,
    get_pdb_tpr,
    copy_mol_to,
    try_copy_to,
)
from locuaz.moleculesutils import get_gro_ziptop_from_pdb, get_gro_ziptop_from_pdb_tleap


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
        target_chains: Iterable,
        binder_chains: Iterable,
        md_config: Dict,
        add_ions: bool = False,
    ) -> "AbstractComplex":
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def from_gro_zip(
        cls,
        *,
        name: str,
        input_dir: Path,
        target_chains: Iterable,
        binder_chains: Iterable,
        gmx_bin: str = "gmx",
    ) -> "AbstractComplex":
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def from_complex(
        cls,
        *,
        name: str,
        branch_path: Path,
        target_chains: Iterable,
        binder_chains: Iterable,
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
    cpt: Optional[FileHandle] = field(kw_only=False, default=None)
    ndx: FileHandle = field(kw_only=False, default=None)

    @classmethod
    def from_pdb(
        cls,
        *,
        name: str,
        input_dir: Union[Path, DirHandle],
        target_chains: Iterable,
        binder_chains: Iterable,
        md_config: Dict,
        add_ions: bool = False,
    ) -> "GROComplex":
        gmx_bin: str = "gmx"
        try:
            in_pdb = Path(input_dir, f"{name}.pdb")
            # This PDB will be backed up and replaced with a fixed up PDB.
            temp_pdb = PDBStructure.from_path(in_pdb)
        except Exception as e:
            print(f"Could not get input PDB file from: {input_dir}", flush=True)
            raise e
        try:
            if md_config["use_tleap"]:
                pdb, gro, top = get_gro_ziptop_from_pdb_tleap(
                    pdb=temp_pdb,
                    target_chains=target_chains,
                    binder_chains=binder_chains,
                )
            else:
                # The PDB should have box info in it, so it'll be pasted onto the GRO. Else,
                # pdb2gmx will compute the binding box of the PDB and use that one.
                pdb, gro, top = get_gro_ziptop_from_pdb(
                    pdb=temp_pdb,
                    target_chains=target_chains,
                    binder_chains=binder_chains,
                    md_config=md_config,
                    add_ions=add_ions
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
            top=top,
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
        target_chains: Iterable,
        binder_chains: Iterable,
        gmx_bin: str = "gmx",
    ) -> "GROComplex":
        try:
            gro = GROStructure.from_path(input_dir / (name + ".gro"))

            top = ZipTopology.from_path_with_chains(
                input_dir / (name + ".zip"),
                target_chains=target_chains,
                binder_chains=binder_chains,
            )
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
            top=top,
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
        branch_path: Path,
        target_chains: Iterable,
        binder_chains: Iterable,
        ignore_cpt: bool = True,
        gmx_bin: str = "gmx",
    ) -> "GROComplex":
        try:
            pdb = PDBStructure.from_path(branch_path / (name + ".pdb"))
            gro = GROStructure.from_path(branch_path / (name + ".gro"))
            top = ZipTopology.from_path_with_chains(
                branch_path / (name + ".zip"),
                target_chains=target_chains,
                binder_chains=binder_chains,
            )

            try:
                traj = XtcTrajectory.from_path(branch_path / (name + ".xtc"))
            except FileNotFoundError:
                try:
                    traj = TrrTrajectory.from_path(branch_path / (name + ".trr"))
                except FileNotFoundError as ee:
                    traj = None
            tpr = TPRFile.from_path(branch_path / (name + ".tpr"))
            if ignore_cpt:
                cpt = None
            else:
                try:
                    cpt = FileHandle(branch_path / (name + ".cpt"))
                except FileNotFoundError as e:
                    cpt = None
        except Exception as ee:
            raise RuntimeError("from_complex() failed.") from ee
        ndx = generate_ndx(
            name,
            pdb=pdb,
            top=top,
        )

        return GROComplex(
            name,
            dir=branch_path,
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
