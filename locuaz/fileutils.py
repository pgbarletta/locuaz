import itertools
import shutil as sh
from functools import singledispatch
from pathlib import Path
from typing import Union

from attrs import define, field


@define
class FileHandle:
    path: Path = field(converter=Path)
    name: str = field(init=False)
    extension: str = field(init=False)

    @path.validator  # type: ignore
    def file_exists(self, attribute, value: Path):
        if not value.is_file():
            raise FileNotFoundError(f"File: {value} doesn't exist.")

    def __attrs_post_init__(self):
        try:
            self.name, self.extension = self.path.name.split(".")
        except ValueError as e:
            self.name = self.path.name
            self.extension = ""
        except Exception as e:
            print(f"Bad input for FileHandle: {self.path}", flush=True)
            raise e

    @classmethod
    def from_existing(cls, name: Path) -> "FileHandle":
        # This method conflicts with the default constructor
        # I'm just leaving it here in case I want to use it later.
        return cls(name)

    def __str__(self) -> str:
        return str(self.path)

    def __fspath__(self) -> str:
        return self.__str__()

    def unlink(self) -> None:
        self.path.unlink()

    def replace_text(self, text_0: str, text_1: str) -> None:
        with open(self.path, "r") as f:
            lineas = f.read()
        lineas = lineas.replace(text_0, text_1)
        with open(self.path, "w") as f:
            f.write(lineas)

    def name_ext(self) -> str:
        return self.path.name


@define
class DirHandle:
    dir_path: Path = field(converter=Path)
    name: str = field(init=False)
    make: bool = field(kw_only=True, default=False)
    force: bool = field(kw_only=True, default=False)
    replace: bool = field(kw_only=True, default=False)

    @dir_path.validator  # type: ignore
    def file_exists(self, attribute, value: Path):
        if not self.make and not value.is_dir():
            raise FileNotFoundError(f"Directory: {value} doesn't exist.")

    def __attrs_post_init__(self):
        self.name = self.dir_path.name
        if self.make:
            try:
                self.dir_path.mkdir()
            except FileExistsError as e_dir_exists:
                if self.force:
                    # Add a numbered prefix to the directory name to avoid conflict.
                    for i in range(1, 100):
                        self.dir_path = Path.joinpath(
                            self.dir_path.parent, str(i) + "-" + self.name
                        )
                        try:
                            Path(self.dir_path).mkdir()
                        except FileExistsError:
                            continue
                        else:
                            self.name = self.dir_path.name
                            break
                    else:
                        print(f"[1:99]-{self.dir_path.name} exist. Can't mkdir.")
                        raise FileExistsError
                elif self.replace:
                    # Delete the conflicting directory.
                    sh.rmtree(self.dir_path)
                    self.dir_path.mkdir()
                    print(f"Replaced dir: {self.dir_path}")
                else:
                    raise e_dir_exists

    def __str__(self) -> str:
        return str(self.dir_path)

    def __fspath__(self) -> str:
        return self.__str__()

    def __truediv__(self, key) -> Union[FileHandle, "DirHandle"]:
        """__truediv__ analog to Path's __truediv__ function, but it also checks the
        existence of the resulting path, whether if its file or dir. Use Path(*args...)
        if you don't want this behaviour

        Raises:
            FileNotFoundError: _description_

        Returns:
            _type_: _description_
        """
        new_path = self.dir_path / key
        if new_path.is_file():
            return FileHandle(new_path)
        elif new_path.is_dir():
            return DirHandle(new_path, make=False)
        else:
            raise FileNotFoundError(f"{new_path} doesn't exist.")


def update_header(file_obj: FileHandle, new_header: str):
    """update_header() overwrites the text file handled changing the first line.

    Args:
        new_header (str): new first line
    """
    with open(file_obj.path, "r") as file:
        texto = file.readlines()
    texto[0] = new_header
    with open(file_obj.path, "w") as file:
        [file.write(linea) for linea in texto]


def catenate(
    out_path: Path, *file_objs: FileHandle, newline_between_files: bool = True
):
    lineas = []
    for file_obj in file_objs:
        with open(file_obj.path, "r") as file:
            lineas.append(file.readlines())
        if newline_between_files:
            lineas.append(["\n"])

    texto = itertools.chain.from_iterable(lineas)
    with open(out_path, "w") as file:
        [file.write(linea) for linea in texto]
    return FileHandle(out_path)


def catenate_pdbs(out_path: Path, *file_objs: FileHandle):
    lineas = []
    for file_obj in file_objs:
        with open(file_obj.path, "r") as file:
            lineas.append(file.readlines())
    texto = itertools.chain.from_iterable(lineas)
    with open(out_path, "w") as file:
        [file.write(linea) for linea in texto if linea[0:3] != "END"]
        file.write("END")
    return FileHandle(out_path)


@singledispatch
def copy_to(obj, dir_path: Union[Path, DirHandle], name=None):
    raise NotImplementedError


@copy_to.register
def _(obj: FileHandle, dir_path: Union[Path, DirHandle], name=None):
    if name is None:
        name = obj.path.name
    new_file = Path(dir_path) / name
    sh.copy(obj.path, new_file)
    return FileHandle(new_file)
