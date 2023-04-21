from typing import Dict, Optional
from pathlib import Path
import shutil as sh

# This will be used to map non-conventional AAs to conventional ones, so the
# scoring functions don't fail.
AA_MAP: Dict[str, str] = {
    "ALA": "ALA",
    "ARG": "ARG",
    "ASN": "ASN",
    "ASP": "ASP",
    "CYS": "CYS",
    "CYX": "CYS",
    "CY2": "CYS",
    "GLN": "GLN",
    "GLU": "GLU",
    "GLY": "GLY",
    "HIS": "HIS",
    "HIE": "HIS",
    "HE1": "HIS",
    "HID": "HIS",
    "ILE": "ILE",
    "LEU": "LEU",
    "LYS": "LYS",
    "MET": "MET",
    "PHE": "PHE",
    "PRO": "PRO",
    "SER": "SER",
    "THR": "THR",
    "TRP": "TRP",
    "TYR": "TYR",
    "VAL": "VAL",
}


class GromacsError(Exception):
    pass

class UserInputError(Exception):
    pass

def launch_biobb(biobb_obj, *, can_write_console_log: bool = False, backup_dict: Optional[Path] = None) -> None:
    """

    :param biobb_obj:
    :type biobb_obj:
    :param can_write_console_log:
    :type can_write_console_log:
    :param backup_dict:
    :type backup_dict:
    """
    biobb_obj.can_write_console_log = can_write_console_log
    err = biobb_obj.launch()
    if err == 0 or err is None:
        try:
            Path(biobb_obj.out_log.name).unlink()
            Path(biobb_obj.err_log.name).unlink()
        except (FileNotFoundError, Exception):
            # log files may get the same name when running parallel biobb processes
            # and one run may try to delete a log file that was already deleted.
            pass
    else:
        log_out = Path(biobb_obj.out_log.name)
        log_err = Path(biobb_obj.err_log.name)
        if backup_dict:
            try:
                new_log_out = backup_dict / f"{log_out.name}_out.txt"
                sh.move(log_out, new_log_out)
                log_out = new_log_out

                new_log_err = backup_dict / f"{log_err.name}_err.txt"
                sh.move(log_err, new_log_err)
                log_err = new_log_err
            except FileNotFoundError:
                # There may be concurrent biobb processes so this file may have been copied already. In this case,
                # the user will probably get a log file in the wrong `backup_dict`, but I can't do much about this
                pass

        raise GromacsError(f"{biobb_obj} failed. Check {log_out} and {log_err} .")


def ext(name: str, suffix: str) -> str:
    """ext utility function, so I don't have to worry about the extension.

    Args:
        name (str): filename with or without the extension.
        suffix (str): desired extension.

    Returns:
        str: filename with the extension.
    """
    return f"{name.split('.')[0]}.{suffix}"
