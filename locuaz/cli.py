import argparse
import glob
import os
import pickle
import shutil as sh
import sys
from warnings import warn, simplefilter
from collections import defaultdict
from itertools import product
from pathlib import Path
from queue import PriorityQueue
from typing import Dict, List, Final, Set, Tuple, Union, Any

import yaml
import networkx as nx

from locuaz.validatore import Validatore
from locuaz.primitives import UserInputError


def get_dir_size(folder: Path) -> float:
    B_TO_MB: Final = 1048576
    total_size = 0
    for path, _, files in os.walk(folder):
        for f in files:
            fp = os.path.join(path, f)
            total_size += os.path.getsize(fp)
    return total_size / B_TO_MB


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

    # Check mode
    if mode != config["main"]["mode"]:
        warn(f"Warning, CLI input {mode} doesn't match {config['main']['mode']}. Overwriting option.")
        config["main"]["mode"] = mode

    config["main"]["debug"] = debug

    # Check input work dir
    if Path(config["paths"]["work"]).is_dir():
        start = False
    else:
        root_dir = Path(config["paths"]["work"]).parent
        assert root_dir.is_dir(), f"Invalid input work dir: {config['paths']['work']}"
        assert (
                mode == "evolve"
        ), "`--mode` is not set to 'evolve', a `work` folder with valid branches is needed."

        start = True

    # Check use of tleap
    if config["md"].get("use_tleap", False):
        assert "tleap" in config["paths"], f"Specify path to tleap files."
        assert (
                get_dir_size(config["paths"]["tleap"]) < 10
        ), f"tleap dir is heavier than 10Mb. Choose a dir with only the necessary tleap files."
        warn(
            "tleap script will be used. Options: 'water_type' and 'force_field' will be ignored.")

    # Check matching lengths of mutating_resSeq and mutating_resname.
    for resnames, resSeqs in zip(config["binder"]["mutating_resname"], config["binder"]["mutating_resSeq"]):
        assert len(resnames) == len(resSeqs), f"mutating_resname({resnames}) and mutating_resSeq ({resSeqs}) have " \
                                              f"different lengths: {len(resnames)} and {len(resSeqs)}, respectively."

    # Check 'branches' and 'generator'
    if ("SPM4" in config["generation"]["generator"]) and (config["protocol"]["new_branches"] > 19):
        raise UserInputError(f"SPM4 mutation generators cannot generate more than 19 branches")

    check_gmxmmpbsa(config)

    # Check consensus threshold and scoring functions
    if config["pruning"]["pruner"] == "consensus":
        assert "consensus_threshold" in config["pruning"], "Set 'consensus_threshold' when using 'consensus' pruner."
        t = config["pruning"]["consensus_threshold"]
        n = len(config["scoring"]["functions"])
        assert t <= n, f"Threshold ({t}) should be equal or lower than the number of scoring functions ({n})"
    elif config["pruning"]["pruner"] == "metropolis":
        assert len(config["scoring"]["functions"]) == 1, f"Can only use 1 scoring function when pruner is 'metropolis'."

    return config, start


def backup_branch(it_fn: Union[str, Path]) -> None:
    it_path = Path(it_fn)
    new_path = Path(it_path.parent, "bu_" + it_path.name)
    warn(f"Found incomplete epoch. Will backup {it_path} to {new_path}")
    sh.move(it_path, new_path)


def append_branches(
        sorted_iters: PriorityQueue, branches: List, prev_epoch: int
) -> str:
    if sorted_iters.empty():
        return ""

    epoch_nbr, branch_str = sorted_iters.get()
    if prev_epoch == 1 or epoch_nbr == prev_epoch:
        branches.append(branch_str)
    else:
        sorted_iters.put((epoch_nbr, branch_str))
        return branch_str

    return append_branches(sorted_iters, branches, epoch_nbr)


def get_valid_branch_dirs(files_and_dirs: List[str], config: Dict) -> List[str]:
    """get_valid_branch_dirs(): filters out paths in input and returns valid branch paths

    Args:
        files_and_dirs (List[str]): list of files and/or dirs
        config (dict): input config.

    Returns:
        List[str]: list of valid branch names
    """

    def is_incomplete(it_path: Path, name: str) -> bool:
        return not Path(it_path, f"{name}.pdb").is_file()

    branch_dirs: List[Path] = []
    incomplete_epochs: Set[str] = set()
    iters_per_epoch: Dict[str, int] = defaultdict(int)
    for filename in files_and_dirs:
        dir_path = Path(config["paths"]["work"], filename)
        if dir_path.is_dir():
            try:
                nbr, *branch_name = Path(dir_path).name.split("-")
            except (Exception,):
                # not an Epoch folder
                continue
            if nbr.isnumeric():
                assert int(nbr) <= config["protocol"]["epochs"], (
                    f"Max of {config['protocol']['epochs']} epochs solicited, but found branch "
                    f"{nbr}-{'-'.join(branch_name)} . Aborting."
                )
                iters_per_epoch[nbr] += 1
                branch_dirs.append(dir_path)
                if is_incomplete(dir_path, config["main"]["name"]):
                    incomplete_epochs.add(nbr)

    valid_iters = []
    for branch_path in branch_dirs:
        nbr, *_ = branch_path.name.split("-")
        if nbr in incomplete_epochs:
            backup_branch(branch_path)
        else:
            valid_iters.append(branch_path.name)

    assert len(valid_iters) > 0, "No valid branches in input."

    return valid_iters


def is_not_epoch_0(branches: Union[List[str], List[Path]], starting_epoch: int):
    for it_fn in branches:
        it_path = Path(it_fn)
        if int(it_path.name.split("-")[0]) == starting_epoch:
            return False
    return True


def lacks_branches(current_branches: Union[List[str], List[Path]],
                   top_branches: Union[List[str], List[Path]], config: Dict[str, Any]) -> bool:
    if config["protocol"]["prevent_fewer_branches"]:
        if config["protocol"]["constant_width"]:
            branches = config["protocol"]["new_branches"]
        else:
            n_top_iters = 1 if len(top_branches) == 0 else len(top_branches)
            branches = config["protocol"]["new_branches"] * n_top_iters
        nbr_branches = len(current_branches)
        starting_epoch = config["main"]["starting_epoch"]

        # Check that the epoch wasn't cut short during its initialization.
        # This may happen if last run was cut during initialize_new_epoch().
        if nbr_branches < branches and is_not_epoch_0(current_branches, starting_epoch):
            for it_fn in current_branches:
                backup_branch(it_fn)
            return True
    return False


def read_key_from_tracking_file(config: Dict[str, Any], tracking: Dict[str, Any], attr: str) -> Dict[str, Any]:
    try:
        config[attr] = tracking[attr]
    except (Exception,):
        warn(f"Could not read {attr} tree from tracking file.")
    return config


def get_tracking_files(config: Dict) -> bool:
    tracking_pkl = Path(config["paths"]["work"], "tracking.pkl")

    try:
        with open(tracking_pkl, "rb") as file:
            tracking: Dict[str, Any] = pickle.load(file)
    except (Exception,):
        warn("Could not read tracking file.")
        return False
    try:
        previous_branches = get_valid_branch_dirs(tracking["previous_branches"], config)
        current_branches = get_valid_branch_dirs(tracking["current_branches"], config)
        top_branches = get_valid_branch_dirs(tracking["top_branches"], config)

        if lacks_branches(current_branches, top_branches, config):
            warn("Will ignore the tracking file. The current branches field is invalid.")
            return False

        config["paths"]["previous_branches"] = previous_branches
        config["paths"]["current_branches"] = current_branches
        config["paths"]["top_branches"] = top_branches
        config["misc"]["epoch_mutated_positions"] = set(tracking["epoch_mutated_positions"])

        if "memory_positions" in config["protocol"]:
            warn(f"Tracking file's 'memory_positions' ignored: {tracking['memory_positions']}.\n"
                 f"Using input config 'memory_positions' instead: {config['protocol']['memory_positions']}")
        else:
            config["protocol"]["memory_positions"] = tracking["memory_positions"]

        if "failed_memory_positions" in config["protocol"]:
            warn(f"Tracking file's 'failed_memory_positions' ignored: {tracking['failed_memory_positions']}.\n"
                 "Using input config 'failed_memory_positions' instead: "
                 f"{config['protocol']['failed_memory_positions']}")
        else:
            config["protocol"]["failed_memory_positions"] = tracking["failed_memory_positions"]

        # Read options related to the DAGs
        for key in ("project_dag", "project_mut_dag", "mutations"):
            config = read_key_from_tracking_file(config, tracking, key)

    except (Exception,):
        warn("Will ignore the tracking file. It is in an invalid state with respect to the working dir.")
        return False

    return True

def set_branches(config: Dict) -> Dict:
    """set_branches Set config["paths"]["current_branches"],
    config["paths"]["previous_branches"] and possibly config["paths"]["top_branches"].

    Args:
        config (Dict): dictionary with input config

    Raises:
        ValueError: when there are no valid branch dirs.
    """
    files_and_dirs = glob.glob(str(Path(config["paths"]["work"], "*")))
    valid_branches = get_valid_branch_dirs(files_and_dirs, config)

    branches: PriorityQueue = PriorityQueue()
    for branch_path in valid_branches:
        nbr, *_ = branch_path.name.split("-")
        branches.put((-int(nbr), branch_path))
    # Branches are sorted by epoch number now
    while not branches.empty():
        config["paths"]["current_branches"] = []
        branch_str = append_branches(branches, config["paths"]["current_branches"], 1)

        if lacks_branches(config["paths"]["current_branches"], [], config):
            # Start over, this time, the incomplete epoch will not be in `iters`.
            continue

        if branch_str != "":
            config["paths"]["previous_branches"] = []
            append_branches(branches, config["paths"]["previous_branches"], 1)
        return config
    raise ValueError("No valid branches in work_dir. Aborting.")


def check_gmxmmpbsa(config: Dict) -> None:
    if config["generation"]["generator"] == "SPM4gmxmmpbsa":
        if "gmxmmpbsa" not in config["scoring"]["functions"]:
            print(f"'gmxmmbpsa' function is necessary to use the 'SPM4mmpbsa' mutator",
                  flush=True, file=sys.stderr)
            raise UserInputError
        else:
            with open(Path(config["paths"]["scoring_functions"], "gmxmmpbsa", "gmxmmpbsa"), 'r') as file:
                for line in file:
                    if "idecomp" in line:
                        return
                    else:
                        continue
                print(f"Input 'gmxmmbpsa' is not performing 'idecomp' residue decomposition. "
                      "This is a prerequisite for the SPM4gmxmmpbsa Mutation Generator.",
                      flush=True, file=sys.stderr)
                raise UserInputError


def set_statistics(config: Dict) -> None:
    """
    set_statistics(): sets a bool in case the user asked for statistics.
    Args:
        config (Dict): input user config

    Returns:
        None
    """
    if "statistics" in config:
        for stat in config["statistics"].keys():
            if stat in {"cmdistance", "interface"}:
                config["statistics"]["get_statistics"] = True
                return
    else:
        config["statistics"] = {"get_statistics": False}


def get_memory(config: Dict) -> Tuple[Set, List[List]]:
    """get_memory() compares character by character of the branch folders resnames to
    find differences among them that would correspond to previously done mutations.
    Small issue in this function: when an epoch was generated from more than 1 top branch,
    the mutations performed on that epoch will be memorized and also the ones that were actually
    performed on a previous epoch and were already present, since it will find differences among them.

    One could remove positions that are repeated in more than 1 memory slot, but it's impossible to know
    if the cause of this was the aforementioned example, or, a previous run with a shorter memory
    (or no memory at all), that mutated the same position more than once in the last N
    (N being equal to `config["protocol"]["memory_size"]`) epochs.
    Hence, no overlap is removed and repetitive resSeqs may be found on the resulting List of Lists.

    Sets: config["misc"]["epoch_mutated_positions"] and config["protocol"]["memory_positions"],


    Args:
        config (Dict): input config options from the user.

    Returns:
        List[List[int]]: resSeqs where different AAs where found. It may be empty, when the branch folders
        being proved don't have the same number of chains and the same length of residues being mutated.
    """

    # Get all the branches sorted by epoch
    iters: PriorityQueue = PriorityQueue()
    for filename in glob.glob(str(Path(config["paths"]["work"], "*"))):
        branch_path = Path(filename)
        if branch_path.is_dir():
            nbr, *_ = Path(branch_path).name.split("-")
            try:
                iters.put((-int(nbr), branch_path))
            except ValueError:
                # not valid branch folder
                continue

    input_n_chains = len(config["binder"]["mutating_chainID"])
    input_n_resSeqs = [len(resSeqs) for resSeqs in config["binder"]["mutating_resSeq"]]
    all_memory_positions: List[List] = []
    # Memorize `config["protocol"]["memory_size"]` epochs.
    for _ in range(config["protocol"]["memory_size"]):
        epoch_iters_paths: List[str] = []
        # Get all the paths for the branches of the current epoch
        append_branches(iters, epoch_iters_paths, 1)

        # Reformat the paths
        epoch_iternames = [
            Path(branch_str).name.split("-")[1:] for branch_str in epoch_iters_paths
        ]
        # Check in case the mutating chains/residues have changed:
        for iterchains in epoch_iternames:
            if len(iterchains) != input_n_chains:
                warn(
                    f"Can't fill requested memory since input 'mutating_chainID' does not "
                    f" match the 'mutating_chainID' previously used on this workspace."
                )
                return set(), [[]]
            old_n_resSeqs = [len(itername[2:]) for itername in iterchains]
            if old_n_resSeqs != input_n_resSeqs:
                warn(
                    f"Can't fill requested memory since input 'mutating_resSeq' does not "
                    f" match the 'mutating_resSeq' previously used on this workspace."
                )
                return set(), [[]]

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
        f'{config["protocol"]["memory_size"]} previous epochs.',
        flush=True,
    )
    # Get the mutated positions on the last epoch
    return set(all_memory_positions[0]), all_memory_positions


def set_empty_dags(config: Dict) -> Dict:
    config["project_dag"] = nx.DiGraph()
    config["project_mut_dag"] = nx.DiGraph()
    return config


def main() -> Tuple[Dict, bool]:
    """Console script for locuaz."""
    simplefilter("default")
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "config_file",
        help="File containing all the necessary parameters to run the protocol",
    )
    parser.add_argument(
        "-m",
        "--mode",
        help="Choose whether to start/restart an evolution protocol or just perform a single task.",
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
    set_statistics(config)
    config["misc"] = {}
    if starts:
        # Set up working dir
        try:
            Path(config["paths"]["work"]).mkdir()
        except FileExistsError as e:
            raise e
    else:
        # Try to read pickle info written by a previous run.
        if get_tracking_files(config):
            print("Read tracking file from work dir.", flush=True)
        else:
            print(
                "Will try to get the previous branches, the current branches and the memory of the "
                "last mutated positions from the work dir. Memory of failed mutations won't be loaded.",
                flush=True)
            config = set_branches(config)
            config = set_empty_dags(config)
            # Set up the memory
            requested_memory = "memory_size" in config["protocol"]
            if requested_memory:
                epoch_mutated_positions, memory_positions = get_memory(config)
                config["misc"]["epoch_mutated_positions"] = epoch_mutated_positions
                has_no_memory = not ("memory_positions" in config["protocol"])
                if has_no_memory:
                    config["protocol"]["memory_positions"] = memory_positions
            else:
                config["misc"]["epoch_mutated_positions"] = set()

    # Set up environment
    os.environ["OMP_PLACES"] = "threads"
    omp_procs_str = str(config["md"]["omp_procs"])
    if os.environ.get("OMP_NUM_THREADS", "") != omp_procs_str:
        warn(f"Setting 'OMP_NUM_THREADS' environment variable to input 'omp_procs' {omp_procs_str}")
        os.environ["OMP_NUM_THREADS"] = omp_procs_str

    return config, starts
