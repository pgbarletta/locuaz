# from biobb_gromacs.gromacs.gmxselect import Gmxselect
# from biobb_analysis.gromacs.gmx_trjconv_str import GMXTrjConvStr
# from biobb_gromacs.gromacs.solvate import Solvate
# from biobb_gromacs.gromacs.grompp import Grompp


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
