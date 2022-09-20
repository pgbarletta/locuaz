import os
import argparse
from typing import Dict, List, Tuple
from pathlib import Path
from queue import PriorityQueue
import glob

import yaml

from validatore import Validatore


def get_raw_config(config_path: str):
    try:
        with open(Path(config_path), "r") as file:
            raw_config = yaml.load(file, yaml.SafeLoader)
    except Exception as e:
        print(
            "Bad .yaml file. Check one of the included configuration files "
            "for a YAML format example.",
            flush=True,
        )
        raise e
    return raw_config


def validate_input(raw_config: Dict, mode: str, debug: bool):
    # Get the cerberus schema located in this same folder
    here_dir = Path(__file__).parent.resolve()
    schema_fn = Path.joinpath(here_dir, "schema.yaml")
    with open(schema_fn, "r") as file:
        schema = yaml.load(file, yaml.SafeLoader)

    # Use the schema to validate the input configuration yaml
    validator = Validatore(schema)
    if not validator.validate(raw_config):
        print(f"Wrong input config file.", flush=True)
        raise ValueError(validator.errors)
    else:
        config = validator.normalized(raw_config)

    if mode != config["main"]["mode"]:
        print(
            f"Warning, CLI input {mode} doesn't match "
            f"{config['main']['mode']}. Overwriting option."
        )
        config["main"]["mode"] = mode
    config["main"]["debug"] = debug

    return config


def append_iterations(
    sorted_iters: PriorityQueue, iterations: List, prev_epoch: int
) -> str:
    if sorted_iters.empty():
        return ""

    epoch_nbr, iter_path = sorted_iters.get()
    if prev_epoch == 1:
        iterations.append(str(iter_path))
    elif epoch_nbr == prev_epoch:
        iterations.append(str(iter_path))
    else:
        return str(iter_path)

    return append_iterations(sorted_iters, iterations, epoch_nbr)


def set_iterations(config: Dict) -> None:
    iters: PriorityQueue = PriorityQueue()
    for filename in glob.glob(str(Path(config["paths"]["work"], "*"))):
        iter_path = Path(filename)
        if iter_path.is_dir():
            nbr, *_ = Path(iter_path).name.split("-")
            iters.put((-int(nbr), iter_path))
    # Iterations are sorted by epoch number now
    config["paths"]["current_iterations"] = []
    iter_str = append_iterations(iters, config["paths"]["current_iterations"], 1)
    if iter_str != "":
        config["paths"]["previous_iterations"] = [iter_str]
        append_iterations(iters, config["paths"]["previous_iterations"], 1)


def main() -> Dict:
    """Console script for locuaz."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "config_file",
        help="File containing all the necessary parameters to run the protocol",
    )
    parser.add_argument(
        "-m",
        "--mode",
        help="Choose wheter to start/restart an evolution protocol or just perform a single task.",
        default="evolve",
        type=str,
        required=False,
        choices=("evolve", "run", "run_npt", "score"),
    )
    parser.add_argument(
        "--debug",
        help="Set/unset debug mode.",
        action="store_true",
    )
    args = parser.parse_args()

    raw_config = get_raw_config(args.config_file)
    config = validate_input(raw_config, args.mode, args.debug)
    if "work" in config["paths"]:
        set_iterations(config)

    # Set up environment
    os.environ["OMP_PLACES"] = "threads"

    return config
