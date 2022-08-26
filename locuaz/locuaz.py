#!/usr/bin/env python

"""Main module."""
import sys
from pathlib import Path
import logging

import cli
import projectutils as pu
from protocol import initialize_new_epoch
from scoring import score, initialize_scoring_folder
from runutils import MDrun
from mututils import MutatorEvoEF2
from run import run_and_score


def main() -> int:

    config = cli.main()
    work_pjct = pu.WorkProject(config)
    this_epoch = work_pjct.epochs[-1]
    iter_name, this_iter = next(iter(this_epoch.items()))

    # mutator = MutatorEvoEF2(work_pjct.config)
    # initialize_new_epoch(work_pjct, mutator)

    if config["main"]["mode"] == "start":
        logging.info(f"Running base NPT run.")
        npt = MDrun.npt(
            this_iter.dir_handle, work_pjct=work_pjct, out_name="npt_" + work_pjct.name
        )
        this_iter.complex = npt(this_iter.complex)

    logging.info(f"Scoring base NPT run.")
    nframes = initialize_scoring_folder(work_pjct, iter_name)
    score(work_pjct, iter_name, nframes)

    if config["main"]["mode"] == "score":
        return 0

    logging.info(f"Starting protocol.")
    mutator = MutatorEvoEF2(work_pjct.config)
    initialize_new_epoch(work_pjct, mutator)
    run_and_score(work_pjct)
    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
