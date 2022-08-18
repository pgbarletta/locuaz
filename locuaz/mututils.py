import typing
from pathlib import Path
from attrs import define, field
from typing import Tuple
from Bio.SeqUtils import seq1, seq3
from projectutils import Epoch, Iteration
from molecules import PDBStructure
from biobb_analysis.gromacs.gmx_trjconv_str import GMXTrjConvStr
from primitives import launch_biobb


@define
class Mutation:

    chainID: str = field(converter=str, kw_only=True)
    resSeq: int = field(converter=int, kw_only=True)
    old_aa: str = field(converter=str, kw_only=True)
    new_aa: str = field(converter=str, kw_only=True)

    def biobb_string(self) -> str:
        """biobb_string returns a string for biobb.model.mutate.mutate

        Returns:
            str: Chain:WT_AA_ThreeLeterCode Resnum MUT_AA_ThreeLeterCode (no whitespace)
        """
        return f"{self.chainID}:{seq3(self.old_aa)}{self.resSeq}{seq3(self.new_aa)}"

    def evoef2_string(self) -> str:
        """evoef2_string returns a string for EvoEF2's --command=BuildMutant

        Returns:
            str: WT_AA_ThreeLeterCode Chain Resnum MUT_AA_ThreeLeterCode (no whitespace)
        """
        return f"{seq3(self.old_aa)}{self.chainID}{self.resSeq}{seq3(self.new_aa)};"


def prepare_old_epoch(epoch: Epoch) -> None:
    """prepare_old_epoch extract 2 PDBs from the npt_`name`.pdb, one with the protein
    and the other with the water and ions. The dry complex will be used for the
    next mutations.

    Args:
        epoch (Epoch): last available epoch, with all its iterations ran, and scored.
    """

    for iter in epoch.values():
        # Protein
        nonwat_pdb = Path(iter.complex.dir_handle) / (
            "nonwat_" + iter.complex.name + ".pdb"
        )

        get_protein = GMXTrjConvStr(
            input_structure_path=str(iter.complex.pdb.file.path),
            input_top_path=str(iter.complex.tpr.file.path),
            output_str_path=str(nonwat_pdb),
            properties={"selection": "Protein"},
        )
        launch_biobb(get_protein)
        iter.nonwat_pdb = PDBStructure.from_path(nonwat_pdb)

        # Water and ions
        wation_pdb = Path(iter.complex.dir_handle) / (
            "wation_" + iter.complex.name + ".pdb"
        )
        get_water_ions = GMXTrjConvStr(
            input_structure_path=str(iter.complex.pdb.file.path),
            input_top_path=str(iter.complex.tpr.file.path),
            output_str_path=str(wation_pdb),
            properties={"selection": "Water_and_ions"},
        )
        launch_biobb(get_water_ions)

        iter.wation_pdb = PDBStructure.from_path(wation_pdb)
        wation_pdb
