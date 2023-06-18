from pathlib import Path
from typing import List
import concurrent.futures as cf
import numpy as np
from numpy.typing import NDArray

from locuaz.projectutils import Branch
from locuaz.basestatistic import BaseStatistic
from locuaz.interface import get_interface_surface


class StatisticInterface(BaseStatistic):
    TIMEOUT_PER_FRAME: int = 2

    def __init__(self, branch: Branch, stats_config: dict) -> None:
        super().__init__(branch, stats_config)

    def __call__(self, start: int, end: int) -> NDArray[float]:
        # The first unused frames will be discarded later
        self.result: NDArray[float] = np.zeros(end)

        if self.nthreads == 1:
            for i in range(start, end):
                _, self.result[i] = get_interface_surface(Path(self.frames_path / f"complex-{i}.pdb"))
        else:
            with cf.ProcessPoolExecutor(max_workers=self.nthreads) as exe:
                futuros: List[cf.Future] = []

                for i in range(start, end):
                    futuros.append(exe.submit(get_interface_surface, Path(self.frames_path / f"complex-{i}.pdb"), i))
                timeout = self.TIMEOUT_PER_FRAME * (end - start)
                try:
                    for futu in cf.as_completed(futuros, timeout=timeout):
                        if futu.exception():
                            print(
                                f"Exception while running {self.name}: {futu.exception()}",
                                flush=True,
                            )
                            raise futu.exception()  # type: ignore
                        j, surface = futu.result()
                        self.result[j] = surface
                except cf.TimeoutError as e:
                    print(f"{self.name} subprocess timed out.", flush=True)
                    raise e

            # Discard the first zeroed frames
            self.result = self.result[start:]
        super().__warn__(str(self), start)
        return self.result

    def __str__(self) -> str:
        return "StatisticInterface"
