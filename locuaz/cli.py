import os
import argparse
from typing import Dict, List, Optional, Set, Tuple, Union, Any
from pathlib import Path
from queue import PriorityQueue
import glob
from itertools import product
from warnings import warn
import shutil as sh
import yaml
import pickle
from collections import defaultdict, deque

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


def backup_iteration(it_fn: Union[str, Path]) -> None:
    it_path = Path(it_fn)
    new_path = Path(it_path.parent, "bu_" + it_path.name)
    # TODO: won't be necessary to cast after 3.9 upgrade
    sh.move(str(it_path), str(new_path))
    warn(f"Found incomplete epoch. Will backup {it_path} to {new_path}")


def append_iterations(
    sorted_iters: PriorityQueue, iterations: List, prev_epoch: int
) -> str:
    if sorted_iters.empty():
        return ""

    epoch_nbr, iter_path = sorted_iters.get()
    if prev_epoch == 1 or epoch_nbr == prev_epoch:
        iterations.append(str(iter_path))
    else:
        sorted_iters.put((epoch_nbr, iter_path))
        return str(iter_path)

    return append_iterations(sorted_iters, iterations, epoch_nbr)


def get_valid_iter_dirs(files_and_dirs: List[str], config: Dict) -> List[Path]:
    """get_valid_iter_dirs Filter out paths in input and return valid iteration paths

    Args:
        files_and_dirs (List[str]): list of files and/or dirs
        name (str): name of the project
        max_epochs (int): max epoch number an iteration may have

    Returns:
        List[Path]: list of valid iteration paths
    """

    def is_incomplete(iter_path: Path, name: str) -> bool:
        return not Path(iter_path, f"{name}.pdb").is_file()

    iter_dirs: List[Path] = []
    incomplete_epochs: Set[str] = set()
    iters_per_epoch: Dict[str, int] = defaultdict(int)
    for filename in files_and_dirs:
        dir_path = Path(filename)
        if dir_path.is_dir():
            try:
                nbr, *iter_name = Path(dir_path).name.split("-")
            except:
                # not an Epoch folder
                continue
            if nbr.isnumeric():
                assert int(nbr) <= config["protocol"]["epochs"], (
                    f"Max of {config['protocol']['epochs']} epochs solicited, but found iteration "
                    f"{nbr}-{'-'.join(iter_name)} . Aborting."
                )
                iters_per_epoch[nbr] += 1
                iter_dirs.append(dir_path)
                if is_incomplete(dir_path, config["main"]["name"]):
                    incomplete_epochs.add(nbr)

    valid_iters = []
    for iter_path in iter_dirs:
        nbr, *_ = iter_path.name.split("-")
        if nbr in incomplete_epochs:
            backup_iteration(iter_path)
        else:
            valid_iters.append(iter_path)

    assert len(valid_iters) > 0, "No valid iterations in input."

    return valid_iters


def is_not_epoch_0(iterations: Union[List[str], List[Path]]):
    for it_fn in iterations:
        it_path = Path(it_fn)
        if it_path.name.split("-")[0] == "0":
            return False
    return True


def lacks_branches(
    current_iterations: Union[List[str], List[Path]],
    branches: int,
    prevent_fewer_branches: bool,
) -> bool:
    if prevent_fewer_branches:
        # Check that the epoch wasn't cut short during its initialization. This may
        # happen if last run was cut during initialize_new_epoch().
        nbr_branches = len(current_iterations)
        if nbr_branches < branches and is_not_epoch_0(current_iterations):
            for it_fn in current_iterations:
                backup_iteration(it_fn)
            return True
    return False


def get_tracking_files(config: Dict) -> bool:
    tracking_pkl = Path(config["paths"]["work"], "tracking.pkl")

    try:
        with open(tracking_pkl, "rb") as file:
            tracking: Dict[str, Any] = pickle.load(file)

        previous_iterations = get_valid_iter_dirs(
            tracking["previous_iterations"], config
        )
        current_iterations = get_valid_iter_dirs(tracking["current_iterations"], config)
        top_iterations = get_valid_iter_dirs(tracking["top_iterations"], config)

        if lacks_branches(
            current_iterations,
            config["protocol"]["branches"],
            config["protocol"]["prevent_fewer_branches"],
        ):
            return False

        config["paths"]["previous_iterations"] = previous_iterations
        config["paths"]["current_iterations"] = current_iterations
        config["paths"]["top_iterations"] = top_iterations
        config["misc"]["epoch_mutated_positions"] = set(
            tracking["epoch_mutated_positions"]
        )
        config["protocol"]["memory_positions"] = deque(tracking["memory_positions"])
        config["protocol"]["failed_memory_positions"] = deque(
            tracking["failed_memory_positions"]
        )

        return True
    except Exception:
        warn(
            "Could not read tracking info. Will try to get the previous iterations, the "
            "current iterations, the memory of the last mutated positions from the work dir. "
            "Memory of failed mutations won't be loaded."
        )
        return False


def set_iterations(config: Dict) -> None:
    """set_iterations Set config["paths"]["current_iterations"],
    config["paths"]["previous_iterations"] and possibly config["paths"]["top_iterations"].

    Args:
        config (Dict): dictionary with input config

    Raises:
        ValueError: when there're no valid iteration dirs.
    """
    files_and_dirs = glob.glob(str(Path(config["paths"]["work"], "*")))
    valid_iters = get_valid_iter_dirs(files_and_dirs, config)

    iters: PriorityQueue = PriorityQueue()
    for iter_path in valid_iters:
        nbr, *_ = iter_path.name.split("-")
        iters.put((-int(nbr), iter_path))
    # Iterations are sorted by epoch number now
    while not iters.empty():
        config["paths"]["current_iterations"] = []
        iter_str = append_iterations(iters, config["paths"]["current_iterations"], 1)

        if lacks_branches(
            config["paths"]["current_iterations"],
            config["protocol"]["branches"],
            config["protocol"]["prevent_fewer_branches"],
        ):
            # Start over, this time, the incomplete epoch will not be in `iters`.
            continue

        if iter_str != "":
            config["paths"]["previous_iterations"] = []
            append_iterations(iters, config["paths"]["previous_iterations"], 1)
        return
    raise ValueError("No valid iterations in work_dir. Aborting.")


def get_memory(config: Dict) -> None:
    """get_memory compare character by character of the iteration folders resnames to
    find differences among them that would correspond to previously done mutations.
    Small issue in this function: when an epoch was generated from more than 1 top iteration,
    the mutations performed on that epoch will be memorized and also the ones that were actually
    performed on a previous epoch and were already present, since it will find differences among them.

    One could remove positions that are repeated in more than 1 memory slot, but it's impossible to know
    if the cause of this was the aforementioned example, or, a previous run with a shorter memory
    (or no memory at all), that mutated the same position more than once in the last N
    (N being equal to `config["protocol"]["memory_size"]`) epochs.
    Hence, no overlap is removedand repetitive resSeqs may be found on the resulting List of Lists.

    Sets: config["misc"]["epoch_mutated_positions"] and config["protocol"]["memory_positions"],


    Args:
        config (Dict): input config options from the user.

    Returns:
        List[List[int]]: resSeqs where different AAs where found. It may be empty, when the iteration folders
        being proved don't have the same number of chains and the same length of residues being mutated.
    """

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
        memory_positions = set()
        for chain_idx in range(input_n_chains):
            resnames = [
                chainID_resname[chain_idx][2:] for chainID_resname in epoch_iternames
            ]

            for str_1, str_2 in list(product(resnames, resnames)):
                for pos, (char_1, char_2) in enumerate(zip(str_1, str_2)):
                    if char_1 != char_2:
                        # Turn the (chain, position) into resSeq
                        memory_positions.add(
                            config["binder"]["mutating_resSeq"][chain_idx][pos]
                        )
        all_memory_positions.append(list(memory_positions))

    print(
        f"Memorized the following positions: {all_memory_positions} from the "
        f'{config["protocol"]["memory_size"]} previous epochs.'
    )
    # Get the mutated positions on the last epoch
    config["misc"]["epoch_mutated_positions"] = set(all_memory_positions[0])
    config["protocol"]["memory_positions"] = all_memory_positions


def define_box_settings(config: Dict) -> Dict:
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
    if starts:
        # Set up working dir
        try:
            Path(config["paths"]["work"]).mkdir()
        except FileExistsError as e:
            raise e
    else:
        # Try to read pickle info written by a previous run.
        if get_tracking_files(config):
            pass
        else:
            set_iterations(config)
            # Set up the memory
            requested_memory = "memory_size" in config["protocol"]
            has_no_memory = not ("memory_positions" in config["protocol"])
            if requested_memory and has_no_memory:
                get_memory(config)
            else:
                config["misc"]["epoch_mutated_positions"] = set()

    config = define_box_settings(config)

    # Set up environment
    os.environ["OMP_PLACES"] = "threads"

    return config, starts
