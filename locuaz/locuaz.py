#!/usr/bin/env python

"""Main module."""
import sys
import logging
from pathlib import Path

import cli
import projectutils as pu
from protocol import initialize_new_epoch
from scoring import score
from run import run_epoch, run_min_nvt_epoch, run_npt_epoch
from prunners import prune


def main() -> int:

    config, start = cli.main()
    log = pu.set_logger(config["main"]["name"], Path(config["paths"]["work"]))
    log.info(f"Setting up work project.")
    work_pjct = pu.WorkProject(config, start)
    log = logging.getLogger(work_pjct.name)
    log.info(f"Launching.")

    if config["main"]["mode"] == "run":
        run_epoch(work_pjct)
        return 0
    elif config["main"]["mode"] == "run_npt":
        run_npt_epoch(work_pjct)
        return 0
    elif config["main"]["mode"] == "score":
        for iter_name, iter in work_pjct.epochs[-1].items():
            log.info(f"Scoring NPT run from iteration: {iter_name}")
            score(work_pjct, iter)
        return 0

    cnt = 0
    while True:
        old_id = work_pjct.epochs[-1].id

        log.info(f"Running epoch {old_id} ({cnt} on this run).")
        run_epoch(work_pjct)

        log.info(f"Scoring epoch {old_id} ({cnt} on this run).")
        for iter_name, iter in work_pjct.epochs[-1].items():
            log.info(f"Scoring NPT run from iteration: {iter_name}")
            score(work_pjct, iter)

        log.info(f"Prunning epoch {old_id}.")
        prune(work_pjct)

        new_id = work_pjct.epochs[-1].id + 1
        cnt += 1
        if new_id > config["protocol"]["epochs"]:
            log.info(f"Done with protocol.")
            break

        log.info(f"Initializing new epoch {new_id} ({cnt} on this run).")
        initialize_new_epoch(work_pjct)

    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
