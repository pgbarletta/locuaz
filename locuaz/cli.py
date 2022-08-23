import argparse
import os
from typing import Dict
import yaml
from pathlib import Path
from cerberus import Validator
import subprocess as sp


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
    config["scoring"]["sf_cnt"] = len(config["scoring"]["functions"])
    if config["scoring"]["sf_cnt"] < config["scoring"]["consensus_threshold"]:
        print(
            f"Wrong input config file. `consensus_threshold` is higher than the number of scoring functions.",
            flush=True,
        )
        raise ValueError

    #
    nchains_chainID = len(config["binder"]["mutating_chainID"])
    nchains_resSeq = len(config["binder"]["mutating_resSeq"])
    if nchains_chainID != nchains_resSeq:
        print(
            f"Wrong input config file. `mutating_chainID` has {nchains_chainID} "
            f"chains, while `mutating_resSeq` has {nchains_resSeq}.  "
            f"They should all have the same length.",
            flush=True,
        )
        raise ValueError

    # If `nprocs` is not set, use as many threads as available
    if config["scoring"]["nprocs"] == 0:
        config["scoring"]["nprocs"] = os.sched_getaffinity(0)
    # If `ngpus` is not set, use as many gpus as available
    if config["md"]["ngpus"] == 0:
        print(
            "Warning: `ngpus` was not set, using `nvidia-smi` to set this "
            "variable, this number may be different that the actual GPUs available "
            "if BinDesigner was launched using a job scheduler.",
            flush=True,
        )
        config["md"]["ngpus"] = int(
            sp.run(
                "nvidia-smi --query-gpu=name --format=csv,noheader | wc -l",
                stdout=sp.PIPE,
                stderr=sp.PIPE,
                shell=True,
                text=True,
            ).stdout.strip()
        )

    # If `branches` is not set, used as many as GPUs are
    if config["main"]["branches"] == 0:
        config["main"]["branches"] = config["md"]["ngpus"]

    # Check that
    necessary_procs = (
        config["md"]["ngpus"] * config["md"]["mpi_procs"] * config["md"]["omp_procs"]
    )
    try:
        available_procs = int(os.getenv("SLURM_CPUS_ON_NODE"))  # type: ignore
    except TypeError as e:
        # TODO: also try to check TORQUE PBS environmental variable.
        available_procs = len(os.sched_getaffinity(0))
    if available_procs < necessary_procs:
        print(
            f"Warning, {config['md']['ngpus']} gpus, {config['md']['mpi_procs']} MPI "
            f"processors and {config['md']['omp_procs']} OMP threads requested. "
            f"{necessary_procs} threads are necessary, but only {available_procs} "
            f"are available.\n Continue only if you know what you're doing.",
            flush=True,
        )

    # This is not strictly necessary yet, but it might eventually be.
    for segment in config["binder"]["mutating_resSeq"]:
        sorted_and_continuous = list(range(segment[0], segment[-1] + 1))
        if sorted_and_continuous != segment:
            print(
                "Wrong input config file. Each list in `mutating_resSeq` "
                "should be a list of contiguous sorted resSeq numbers. Eg: "
                "[ [86, 87, 88], [90, 91, 92] ]",
                flush=True,
            )
            raise ValueError

    return config


def assert_dir(dir: str):
    dir_path = Path(dir)
    assert dir_path.is_dir(), f"Error on input. {dir_path} is not an existing folder."
    return dir_path


def set_start_or_restart(config: Dict) -> Dict:

    # check that the folders actually exist
    for name, input_str in config["paths"].items():
        if name == "input":
            # This one is a list
            for each_input_str in input_str:
                assert_dir(each_input_str)
        else:
            assert_dir(input_str)
    if config["main"]["mode"] != "start":
        epoch_nbrs = []
        for iteration_str in config["paths"]["input"]:
            try:
                # I'm assuming the iteration folder's name wasn't changed and starts
                # as `i-...`, where `i` is the number of the last ran iteration.
                epoch_nbrs.append(int(Path(iteration_str).name.split("-")[0]))
            except ValueError as e:
                print(
                    f"{iteration_str} is not a valid starting iteration folder. ",
                    flush=True,
                )
                raise e
        config["epoch_nbr"] = max(epoch_nbrs)
    config["md"]["gmx"] = Path.joinpath(Path(config["paths"]["gmxrc"]), "gmx")

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
            "Bad .yaml file. Check one of the included configuration files "
            "for a YAML format example.",
            flush=True,
        )
        raise e

    config = validate_config(raw_config)
    config = set_start_or_restart(config)

    return config
