import os
import argparse
from typing import Dict
from pathlib import Path

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

    # Set up environment
    os.environ["OMP_PLACES"] = "threads"

    return config
