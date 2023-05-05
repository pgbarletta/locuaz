import logging
import sys
from pathlib import Path

import locuaz.cli
from locuaz.projectutils import set_logger, WorkProject
from locuaz.epochinitializer import initialize_new_epoch
from locuaz.pruning import prune
from locuaz.run import run_epoch
from locuaz.scoring import score


def main() -> int:
    config, start = locuaz.cli.main()
    log = set_logger(config["main"]["name"], Path(config["paths"]["work"]))
    log.info("-----------------------------")
    log.info("Setting up work project.")
    work_pjct = WorkProject(config, start)
    log = logging.getLogger(work_pjct.name)
    log.info("Launching.")

    if config["main"]["mode"] == "run":
        run_epoch(work_pjct)
        return 0
    elif config["main"]["mode"] == "score":
        try:
            prev_id = work_pjct.epochs[-2].id
            prev_epoch = work_pjct.epochs[-2]
            for iter_name, iteration in prev_epoch.items():
                log.info(f"Scoring NPT run from iteration: {prev_id}-{iter_name}")
                score(work_pjct, iteration)
        except (Exception,):
            log.info("Cannot score previous epoch.")
            pass

        curr_id = work_pjct.epochs[-1].id
        for iter_name, iteration in work_pjct.epochs[-1].items():
            log.info(f"Scoring NPT run from iteration: {curr_id}-{iter_name}")
            score(work_pjct, iteration)
        return 0

    cnt = 0
    while True:
        old_id = work_pjct.epochs[-1].id

        log.info(f"Running epoch {old_id} ({cnt} on this run).")
        run_epoch(work_pjct)

        log.info(f"Scoring epoch {old_id} ({cnt} on this run).")
        for iter_name, iteration in work_pjct.epochs[-1].items():
            log.info(f"Scoring NPT run from iteration: {iter_name}")
            score(work_pjct, iteration)

        log.info(f"Pruning epoch {old_id}.")
        prune(work_pjct)

        new_id = work_pjct.epochs[-1].id + 1
        cnt += 1
        if new_id > config["protocol"]["epochs"]:
            log.info("Done with protocol.")
            break

        log.info(f"Initializing new epoch {new_id} ({cnt} on this run).")
        initialize_new_epoch(work_pjct, log)
    print("Done with protocol")
    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
