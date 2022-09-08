import argparse
import logging
from typing import Dict
from pathlib import Path
from functools import singledispatch
import subprocess as sp
import os
from collections.abc import Sequence

from cerberus import Validator
import yaml

from fileutils import FileHandle

# TODO: move validate_config() to this class
class Valida(Validator):
    def _validate_contains_any_of(self, constraints, field, value):
        """_validate_contains_any_of

        The rule's arguments are validated against this schema:
        {'type': 'list'}
        """
        if constraints:
            if len(set(value.keys()).intersection(constraints)) == 0:
                self._error(field, f"Must contain any of: {constraints}")

    def _validate_step_bigger_than(self, other, field, value):
        """_validate_step_bigger_than

        The rule's arguments are validated against this schema:
        {'type': 'string'}
        """
        if other:
            n = self.document[other]
            if n > value[1]:
                self._error(
                    field,
                    f"step ({value[1]}) is lower than `{other}`({n}). "
                    "This would dedicate threads to more than 1 run.",
                )

    def _validate_sorted(self, flag, field, value):
        """_validate_sorted

        The rule's arguments are validated against this schema:
        {'type': 'boolean'}
        """
        if flag:
            if sorted(value) != value:
                self._error(field, "should be an incrementally sorted list.")


def validate_config(config: Dict) -> Dict:

    # Get the cerberus schema located in this same folder
    here_dir = Path(__file__).parent.resolve()
    schema_fn = Path.joinpath(here_dir, "schema.yaml")
    with open(schema_fn, "r") as file:
        schema = yaml.load(file, yaml.SafeLoader)

    # Use the schema to validate the input configuration yaml
    validator = Valida(schema)
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

    if config["scoring"]["nprocs"] > 64:
        print(
            f'Warning: {config["scoring"]["nprocs"]} threads requested for '
            "scoring. Make sure you have enough RAM for running as many parallel jobs. ",
            flush=True,
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

    if ("SPM" in config["main"]["mutator"]) and (config["main"]["branches"] > 19):
        print(
            f"{config['main']['branches']} is over 19 but this mutator generates "
            "mutations for 1 position, hence, 19 mutations is the maximum number "
            "of possible mutations. Aborting.",
            flush=True,
        )
        raise ValueError

    return config


@singledispatch
def assert_dir(dir) -> Path:
    raise NotImplementedError


@assert_dir.register
def _(dir: str) -> Path:
    dir_path = Path(dir)
    assert dir_path.is_dir(), f"Error on input. {dir_path} is not an existing folder."
    return dir_path


@assert_dir.register
def _(dir: Sequence) -> Path:
    for each_dir_str in dir:
        dir_path = assert_dir(each_dir_str)
    return dir_path


def check_input_paths(config: Dict) -> Dict:

    # check that the folders actually exist
    for input_str in config["paths"].values():
        assert_dir(input_str)
    if "root" in config["paths"]:
        print(f"Starting in {config['paths']['root']}", flush=True)
    else:
        if "previous_iterations" not in config["paths"]:
            config["paths"]["previous_iterations"] = []
        if "current_iterations" not in config["paths"]:
            raise ValueError(
                "`current_iterations` must be specified " "if there is no `root` path."
            )
        print(f"Restarting from {config['paths']['current_iterations']}", flush=True)
        for iteration_str in (
            config["paths"]["previous_iterations"]
            + config["paths"]["current_iterations"]
        ):
            nbr, *chains_resnames = Path(iteration_str).name.split("-")
            if not nbr.isnumeric():
                logging.error(
                    f"{iteration_str} is not a valid iteration folder."
                    f"{nbr} is not a valid epoch number"
                )
                raise ValueError
            for chain_resname in chains_resnames:
                chainID, resname = chain_resname.split("_")
                if not (len(chainID) == 1 and chainID.isupper()):
                    logging.error(
                        f"{iteration_str} is not a valid iteration folder."
                        f"{chainID} is not a valid chainID number"
                    )
                    raise ValueError
                # Don't have checks for resname yet.
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
        help="Wheter to start/restart a protocol or just perform a single task.",
        default="",
        type=str,
        required=False,
        choices=["protocol", "run", "run_npt", "score"],
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
    config = check_input_paths(config)

    if args.mode != "" and args.mode != config["main"]["mode"]:
        print(
            f"Warning, CLI input {args.mode} doesn't match "
            f"{config['main']['mode']}. Overwriting option."
        )
        config["main"]["mode"] = args.mode

    # Set up environment
    os.environ["OMP_PLACES"] = "threads"

    return config
