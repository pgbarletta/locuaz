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

    if not work_pjct.epochs[-1].nvt_done:
        run_min_nvt_epoch(work_pjct)
        work_pjct.epochs[-1].nvt_done = True

    if not work_pjct.epochs[-1].npt_done:
        run_npt_epoch(work_pjct)
        work_pjct.epochs[-1].npt_done = True

    for cnt in range(config["protocol"]["epochs"] + 1):
        old_id = work_pjct.epochs[-1].id

        log.info(f"Scoring epoch {old_id} ({cnt} on this run).")
        for iter_name, iter in work_pjct.epochs[-1].items():
            log.info(f"Scoring NPT run from iteration: {iter_name}")
            score(work_pjct, iter)

        log.info(f"Prunning epoch {old_id}.")
        prune(work_pjct)

        top_itrs_str = " ; ".join(
            [
                f"{iter.epoch_id}-{iter.iter_name}"
                for iter in work_pjct.epochs[-1].top_iterations.values()
            ]
        )
        log.info(f"Top iterations: {top_itrs_str}")

        if cnt == config["protocol"]["epochs"]:
            log.info(f"Done with protocol.")
            break

        log.info(f"Initializing new epoch {old_id+1} ({cnt+1} on this run).")
        initialize_new_epoch(work_pjct)

        log.info(f"Running epoch {old_id+1} ({cnt+1} on this run).")
        run_epoch(work_pjct)

    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
