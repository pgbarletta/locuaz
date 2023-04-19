import logging
from typing import Dict, Optional, Type

from .projectutils import Iteration
from .basestatistic import BaseStatistic
from .statisticcmdistance import StatisticCMDistance
from .statisticinterface import StatisticInterface

all_stats: Dict[str, Type[BaseStatistic]] = {
    "cmdistance": StatisticCMDistance,
    "interface": StatisticInterface,
}


def run_stats(iteration: Iteration, config: dict, start: int, end: int,
              log: Optional[logging.Logger] = None) -> None:
    if config["statistics"]["get_statistics"]:
        for stat, stat_opts in config["statistics"].items():
            if stat == "get_statistics":
                continue
            log.info(f"Getting stat '{stat}' from the trajectory. ", extra={'end': ''})
            run_stat = all_stats[stat](iteration, stat_opts)
            stat_result = run_stat(start, end)
            iteration.set_stat(stat, stat_result, log)

        iteration.write_down_stats()