import os
import argparse
from typing import Dict, List, Optional, Set, Tuple
from pathlib import Path
from queue import PriorityQueue
import glob
from itertools import product
from warnings import warn
import shutil as sh

import yaml

from validatore import Validatore


def get_raw_config(config_path: str):
    try:
        with open(Path(config_path), "r") as file:
            raw_config = yaml.load(file, yaml.SafeLoader)
    except Exception as e:
        raise ValueError(
            "Bad config YAML file. Check one of the example config files."
        ) from e
    return raw_config


def validate_input(raw_config: Dict, mode: str, debug: bool) -> Tuple[Dict, bool]:
    # Get the cerberus schema located in this same folder
    here_dir = Path(__file__).parent.resolve()
    schema_fn = Path.joinpath(here_dir, "schema.yaml")
    with open(schema_fn, "r") as file:
        schema = yaml.load(file, yaml.SafeLoader)

    # Use the schema to validate the input configuration yaml
    validator = Validatore(schema)
    if not validator.validate(raw_config):  # type: ignore
        raise ValueError(validator.errors)  # type: ignore
    else:
        config = validator.normalized(raw_config)  # type: ignore

    if mode != config["main"]["mode"]:
        warn(
            f"Warning, CLI input {mode} doesn't match {config['main']['mode']}."
            f"Overwriting option."
        )
        config["main"]["mode"] = mode

    config["main"]["debug"] = debug

    if Path(config["paths"]["work"]).is_dir():
        start = False
    else:
        root_dir = Path(config["paths"]["work"]).parent
        assert root_dir.is_dir(), f"Invalid input work dir: {config['paths']['work']}"
        assert (
            mode == "evolve"
        ), "`--mode` is not set to 'evolve', a `work` folder with vaild iterations is needed."

        start = True

    return config, start


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
        sorted_iters.put((epoch_nbr, iter_path))
        return str(iter_path)

    return append_iterations(sorted_iters, iterations, epoch_nbr)


def is_incomplete(iter_path: Path, name: str) -> bool:
    return not Path(iter_path, f"{name}.pdb").is_file()


def list_iteration_dirs(wrk_path: Path, name: str) -> List[Path]:
    iter_dirs: List[Path] = []
    incomplete_epochs: Set[str] = set()
    for filename in glob.glob(str(Path(wrk_path, "*"))):
        iter_path = Path(filename)
        if iter_path.is_dir():
            try:
                nbr, *_ = Path(iter_path).name.split("-")
                if nbr.isnumeric():
                    iter_dirs.append(iter_path)
                    if is_incomplete(iter_path, name):
                        incomplete_epochs.add(nbr)
            except:
                # not an Epoch folder
                pass

    valid_iters = []
    for iter_path in iter_dirs:
        nbr, *_ = iter_path.name.split("-")
        if nbr in incomplete_epochs:
            new_path = Path(iter_path.parent, "bu_" + iter_path.name)
            # TODO: won't be necessary to cast after 3.9 upgrade
            sh.move(str(iter_path), str(new_path))
            warn(f"Found incomplete epoch. Will backup {iter_path} to {new_path}")
        else:
            valid_iters.append(iter_path)
    return valid_iters


def set_iterations(config: Dict) -> None:

    valid_iters = list_iteration_dirs(
        Path(config["paths"]["work"]), config["main"]["name"]
    )

    iters: PriorityQueue = PriorityQueue()
    for iter_path in valid_iters:
        nbr, *_ = iter_path.name.split("-")
        iters.put((-int(nbr), iter_path))

    # Iterations are sorted by epoch number now
    config["paths"]["current_iterations"] = []
    iter_str = append_iterations(iters, config["paths"]["current_iterations"], 1)
    if iter_str != "":
        config["paths"]["previous_iterations"] = []
        append_iterations(iters, config["paths"]["previous_iterations"], 1)


def remove_overlap(v: List):
    not_overlapped = []
    for i, j in zip(v[:-1], v[1:]):
        not_overlapped.append(list(set(i) - set(j)))
    return not_overlapped + [v[-1]]


# Small bug in this function: when an epoch was generated from more than 1 top iteration
# The mutations performed on that epoch will be memorized and also the ones that were actually
# performed on a previous epoch, since it will find differences among them
def get_memory(config: Dict) -> Optional[List[List[int]]]:

    # Get all the iterations sorted by epoch
    iters: PriorityQueue = PriorityQueue()
    for filename in glob.glob(str(Path(config["paths"]["work"], "*"))):
        iter_path = Path(filename)
        if iter_path.is_dir():
            nbr, *_ = Path(iter_path).name.split("-")
            try:
                iters.put((-int(nbr), iter_path))
            except ValueError:
                # not valid iteration folder
                continue

    input_n_chains = len(config["binder"]["mutating_chainID"])
    input_n_resSeqs = [len(resSeqs) for resSeqs in config["binder"]["mutating_resSeq"]]
    all_memory_positions = []
    # Memorize `config["protocol"]["memory_size"]` epochs.
    for _ in range(config["protocol"]["memory_size"]):
        epoch_iters_paths: List[str] = []
        # Get all the paths for the iterations of the current epoch
        append_iterations(iters, epoch_iters_paths, 1)

        # Reformat the paths
        epoch_iternames = [
            Path(iter_str).name.split("-")[1:] for iter_str in epoch_iters_paths
        ]
        # Check in case the mutating chains/residues have changed:
        for iterchains in epoch_iternames:
            if len(iterchains) != input_n_chains:
                warn(
                    f"Can't fill requested memory since input 'mutating_chainID' does not "
                    f" match the 'mutating_chainID' previously used on this workspace."
                )
                return None
            old_n_resSeqs = [len(itername[2:]) for itername in iterchains]
            if old_n_resSeqs != input_n_resSeqs:
                warn(
                    f"Can't fill requested memory since input 'mutating_resSeq' does not "
                    f" match the 'mutating_resSeq' previously used on this workspace."
                )
                return None

        # Find the positions and chains that were mutated in this epoch
        different_positions = set()
        for i in range(input_n_chains):
            resname = [chainID_resname[i][2:] for chainID_resname in epoch_iternames]

            for str_1, str_2 in list(product(resname, resname)):
                for j, (char_1, char_2) in enumerate(zip(str_1, str_2)):
                    if char_1 != char_2:
                        different_positions.add((i, j))

        # Turn the (chain, position) into resSeq
        memory_positions = []
        for chain_idx, pos in different_positions:
            memory_positions.append(config["binder"]["mutating_resSeq"][chain_idx][pos])
        all_memory_positions.append(memory_positions)

    fix_all_memory_positions = remove_overlap(all_memory_positions)
    print(
        f"Memorized the following positions: {fix_all_memory_positions} from the "
        f'{config["protocol"]["memory_size"]} previous epochs.'
    )

    return fix_all_memory_positions


def misc_settings(config: Dict) -> Dict:
    # config["md"]["box"] overrides config["md"]["dist_to_box"]
    if "box" in config["md"]:
        config["md"]["dist_to_box"] = None
        config["md"]["box"] = " ".join([str(length) for length in config["md"]["box"]])
    else:
        if "dist_to_box" in config["md"]:
            pass
        else:
            config["md"]["dist_to_box"] = 1.0
    return config


def main() -> Tuple[Dict, bool]:
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
    config, starts = validate_input(raw_config, args.mode, args.debug)
    set_iterations(config)

    requested_memory = "memory_size" in config["protocol"]
    has_no_memory = not ("memory_positions" in config["protocol"])
    if not starts and requested_memory and has_no_memory:
        memory_positions = get_memory(config)
        config["protocol"]["memory_positions"] = memory_positions
    config = misc_settings(config)

    # Set up environment
    os.environ["OMP_PLACES"] = "threads"

    return config, starts
