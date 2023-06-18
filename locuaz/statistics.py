import logging
from typing import Dict, Optional, Type

from locuaz.projectutils import Branch
from locuaz.basestatistic import BaseStatistic
from locuaz.statisticcmdistance import StatisticCMDistance
from locuaz.statisticinterface import StatisticInterface

all_stats: Dict[str, Type[BaseStatistic]] = {
    "cmdistance": StatisticCMDistance,
    "interface": StatisticInterface,
}


def run_stats(branch: Branch, config: dict, start: int, end: int,
              log: Optional[logging.Logger] = None) -> None:
    if config["statistics"]["get_statistics"]:
        for stat, stat_opts in config["statistics"].items():
            if stat == "get_statistics":
                continue
            log.info(f"Getting stat '{stat}' from the trajectory. ", extra={'end': ''})
            run_stat = all_stats[stat](branch, stat_opts)
            stat_result = run_stat(start, end)
            branch.set_stat(stat, stat_result, log)

        branch.write_down_stats()
