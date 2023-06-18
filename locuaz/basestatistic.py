import warnings
from pathlib import Path
from typing import Dict, Tuple, Callable
from abc import abstractmethod

import numpy as np
from numpy.typing import NDArray

from locuaz.projectutils import Branch


class BaseStatistic:
    frames_path: Path
    trj_suffix: str
    name: str
    result: NDArray[float]
    available_warnings: Tuple[str]
    visitor: Dict[str, Callable]
    stats_config: Dict[str, float]

    def __init__(self, branch: Branch, stats_config: Dict) -> None:
        self.frames_path = Path(branch.score_dir)
        self.trj_suffix = Path(branch.complex.tra).suffix
        self.name = branch.complex.name
        self.nthreads = stats_config.get("nthreads", 1)

        self.visitor = {"warn_above": self.__warn_above__,
                        "warn_below": self.__warn_below__,
                        "warn_above_relative": self.__warn_above_relative__,
                        "warn_below_relative": self.__warn_below_relative__,
                        "warn_variance": self.__warn_variance__}

        for miss_opt in [ opt for opt in self.visitor.keys() if opt not in stats_config ]:
            del self.visitor[miss_opt]
        self.stats_config = stats_config

    @abstractmethod
    def __call__(self, start: int, end: int) -> NDArray[float]:
        pass

    def __warn__(self, stat_name: str, offset: int) -> None:
        assert len(self.result) > 0, f"ERROR. {stat_name} failed to get any results."
        for opt, func in self.visitor.items():
            threshold: float = self.stats_config[opt]
            warn, frames = func(threshold, offset)
            if warn:
                if len(frames) > 0:
                    warnings.warn(
                        f"{stat_name}: {opt} found frames above the threshold {threshold}: {frames}")
                else:
                    warnings.warn(
                        f"{stat_name}: {opt} found the trajectory to be above the threshold {threshold}.")

    def __warn_above__(self, threshold: float, offset: int) -> Tuple[bool, NDArray]:
        failing_frames = self.result > threshold
        return sum(failing_frames) > 0, list(np.where(failing_frames)[0] + offset)

    def __warn_below__(self, threshold: float, offset: int) -> Tuple[bool, NDArray]:
        failing_frames = self.result < threshold
        return sum(failing_frames) > 0, list(np.where(failing_frames)[0] + offset)

    def __warn_above_relative__(self, threshold: float, offset: int) -> Tuple[bool, NDArray]:
        first_val = threshold * self.result[0]
        failing_frames = self.result > first_val
        return sum(failing_frames) > 0, list(np.where(failing_frames)[0] + offset)

    def __warn_below_relative__(self, threshold: float, offset: int) -> Tuple[bool, NDArray]:
        first_val = threshold * self.result[0]
        failing_frames = self.result < first_val
        return sum(failing_frames) > 0, list(np.where(failing_frames)[0] + offset)

    def __warn_variance__(self, threshold: float, offset: int) -> Tuple[bool, NDArray]:
        return np.std(self.result) > threshold, []
