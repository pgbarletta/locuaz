from typing import Dict

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


def launch_biobb(biobb_obj, can_write_console_log: bool = False) -> None:
    biobb_obj.can_write_console_log = can_write_console_log
    err = biobb_obj.launch()
    assert (
        err is 0 or err is None
    ), f"{biobb_obj} failed. Check {biobb_obj.out_log} and {biobb_obj.err_log}."


def ext(name: str, ext: str) -> str:
    """ext utility function so I don't have to worry about the extension.

    Args:
        name (str): filename with or without the extension.

    Returns:
        str: filename with the extension.
    """
    return f"{name.split('.')[0]}.{ext}"
