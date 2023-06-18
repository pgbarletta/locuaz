from typing import Dict
import numpy as np

from numpy.typing import NDArray
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array

from locuaz.projectutils import Branch
from locuaz.molecules import XtcTrajectory, PDBStructure
from locuaz.basestatistic import BaseStatistic


class StatisticCMDistance(BaseStatistic):
    pdb_top: PDBStructure
    trj: XtcTrajectory

    def __init__(self, branch: Branch, stats_config: Dict) -> None:
        super().__init__(branch, stats_config)
        self.pdb_top = XtcTrajectory.from_path(self.frames_path / f"fix_{self.name}.pdb")
        self.trj = XtcTrajectory.from_path(self.frames_path / f"fix_{self.name}{self.trj_suffix}")

    def __call__(self, start: int, end: int) -> NDArray[float]:
        u = mda.Universe(str(self.pdb_top), str(self.trj))

        # Assumes the convention target=A, binder=B
        target = u.select_atoms("chainID A or segid A")
        binder = u.select_atoms("chainID B or segid B")

        # Calculate the distance between the centers of mass
        self.result = np.array([distance_array(target.center_of_mass(), binder.center_of_mass()).item()
                                for _ in u.trajectory[start:end]])

        super().__warn__(str(self), start)
        return self.result

    def __str__(self) -> str:
        return "StatisticCenterOfMassDistance"
