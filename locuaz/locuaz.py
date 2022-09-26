#!/usr/bin/env python

"""Main module."""
import sys
import logging

import cli
import projectutils as pu
from protocol import initialize_new_epoch
from scoring import score
from run import run_epoch, run_min_nvt_epoch, run_npt_epoch
from prunners import prune


def main() -> int:

    config = cli.main()
    work_pjct = pu.WorkProject(config)
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

    if not work_pjct.epochs[-1].npt_started:
        run_min_nvt_epoch(work_pjct)
        work_pjct.epochs[-1].npt_started = True

    if not work_pjct.epochs[-1].npt_finished:
        run_npt_epoch(work_pjct)
        work_pjct.epochs[-1].npt_finished = True

    for cnt in range(config["protocol"]["epochs"]):
        epoch_id = work_pjct.epochs[-1].id

        log.info(f"Scoring epoch {epoch_id} ({cnt} on this run).")
        for iter_name, iter in work_pjct.epochs[-1].items():
            log.info(f"Scoring NPT run from iteration: {iter_name}")
            score(work_pjct, iter)

        log.info(f"Prunning epoch {epoch_id}.")
        prune(work_pjct)
        log.info(f"Top iterations: {work_pjct.epochs[-1].top_iterations.keys()}.")

        log.info(f"Initializing new epoch {epoch_id+1} ({cnt+1} on this run).")
        initialize_new_epoch(work_pjct)

        log.info(f"Running epoch {epoch_id} ({cnt} on this run).")
        run_epoch(work_pjct)

    log.info(f"Done with protocol.")

    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
