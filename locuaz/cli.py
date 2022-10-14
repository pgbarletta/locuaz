import os
import argparse
from typing import Dict, List, Optional
from pathlib import Path
from queue import PriorityQueue
import glob
from itertools import product


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
    if not validator.validate(raw_config):  # type: ignore
        print(f"Wrong input config file.", flush=True)
        raise ValueError(validator.errors)  # type: ignore
    else:
        config = validator.normalized(raw_config)  # type: ignore

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
        sorted_iters.put((epoch_nbr, iter_path))
        return str(iter_path)

    return append_iterations(sorted_iters, iterations, epoch_nbr)


def set_iterations(config: Dict) -> None:
    if "work" not in config["paths"]:
        return
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
        config["paths"]["previous_iterations"] = []
        append_iterations(iters, config["paths"]["previous_iterations"], 1)


def remove_overlap(v: List):
    not_overlapped = []
    for i, j in zip(v[:-1], v[1:]):
        not_overlapped.append(list(set(i) - set(j)))
    return not_overlapped + [v[-1]]


# Small bug in this function: when an epoch was generated from more than 1 top iteration
# The mutations performed on that epoch will be memorized and also the ones that were actually
# performed on a previous epoch, since it will find differences among
def get_memory(config: Dict) -> Optional[List[List[int]]]:
    starts = "root" in config["paths"]
    requested_memory = config["protocol"].get("memory_size", False)
    has_set_memory = config["protocol"].get("memory_positions", False)
    if starts or not requested_memory or has_set_memory:
        return None

    # Get all the iterations sorted by epoch
    iters: PriorityQueue = PriorityQueue()
    for filename in glob.glob(str(Path(config["paths"]["work"], "*"))):
        iter_path = Path(filename)
        if iter_path.is_dir():
            nbr, *_ = Path(iter_path).name.split("-")
            iters.put((-int(nbr), iter_path))

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
                print(
                    f"Can't fill requested memory since input 'mutating_chainID' does not "
                    f" match the 'mutating_chainID' previously used on this workspace."
                )
                return None
            old_n_resSeqs = [len(itername[2:]) for itername in iterchains]
            if old_n_resSeqs != input_n_resSeqs:
                print(
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
    set_iterations(config)
    memory_positions = get_memory(config)
    if memory_positions:
        config["protocol"]["memory_positions"] = memory_positions
    config = misc_settings(config)

    # Set up environment
    os.environ["OMP_PLACES"] = "threads"

    return config
