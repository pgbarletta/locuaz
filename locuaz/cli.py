import argparse
import sys
from typing import Dict
import yaml
from pathlib import Path
from enum import Enum
from cerberus import Validator  # type: ignore
from fileutils import DirHandle  # type: ignore


def validate_config(config: Dict) -> Dict:

    # Get the cerberus schema located in this same folder
    here_dir = Path(__file__).parent.resolve()
    schema_fn = Path.joinpath(here_dir, "schema.yaml")
    with open(schema_fn, "r") as file:
        schema = yaml.load(file, yaml.SafeLoader)

    # Use the schema to validate the input configuration yaml
    validator = Validator(schema)
    if not validator.validate(config):
        print(f"Wrong input config file.", flush=True)
        raise ValueError(validator.errors)
    else:
        config = validator.normalized(config)

    # Perform extra checks on the input config file.
    config["scoring"]["sf_cnt"] = len(config["scoring"]["scoring_functions"])
    if config["scoring"]["sf_cnt"] < config["scoring"]["consensus_threshold"]:
        print(
            f"Wrong input config file. `consensus_threshold` is higher than the number of scoring functions.",
            flush=True,
        )
        raise ValueError

    for seq_a, seq_b in zip(
        config["binder"]["peptide_reference"], config["binder"]["residues_mod"]
    ):
        if len(seq_a) != len(seq_a):
            print(
                f"Wrong input config file. `peptide_reference` and `residues_mod` should have the same length.",
                flush=True,
            )
            raise ValueError

    return config


def assert_dir(dir: str):
    dir_path = Path(dir)
    assert dir_path.is_dir(), f"Error on input. {dir_path} is not an existing folder."
    return dir_path


def set_start_or_restart(config: Dict) -> Dict:

    if config["paths"]["iteration"] == "start":
        for name, input_str in config["paths"].items():
            if name == "iteration":
                continue
            assert_dir(input_str)
    else:
        dir_path = assert_dir(config["paths"]["iteration"])
        try:
            # I'm assuming the iteration folder's name wasn't changed and starts
            # as `i-...`, where `i` is the number of the last ran iteration.
            config["epoch_nbr"] = int(dir_path.name.split("-")[0])
            config["paths"]["data"] = config["paths"]["iteration"]
            config["paths"]["input"] = config["paths"]["iteration"]
        except ValueError as e:
            print(f"{dir_path} is not a valid starting iteration folder. ", flush=True)
            raise e

    return config


def main() -> Dict:
    """Console script for locuaz."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "config_file",
        help="File containing all the necessary parameters to run the protocol",
    )
    args = parser.parse_args()

    try:
        with open(Path(args.config_file), "r") as file:
            raw_config = yaml.load(file, yaml.SafeLoader)
    except Exception as e:
        print(
            "Bad .yaml file. Check one of the included configuration files for a YAML format example.",
            flush=True,
        )
        raise e

    config = validate_config(raw_config)
    config = set_start_or_restart(config)

    return config
