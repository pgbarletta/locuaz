from pathlib import Path
from attrs import define, field, validators
from typing import (
    Iterable,
    Iterator,
    List,
    Sequence,
    Set,
    Dict,
    Tuple,
    Union,
    Deque,
    Optional,
)
import logging
from statistics import mean, stdev
from collections import deque
from collections.abc import MutableMapping
from queue import PriorityQueue
import shutil as sh
import time
import pickle
from warnings import warn

from Bio.SeqUtils import seq1
import numpy as np
import MDAnalysis as mda

from fileutils import FileHandle, DirHandle, copy_to
from complex import AbstractComplex, GROComplex
from abstractscoringfunction import AbstractScoringFunction
from scoringfunctions import scoringfunctions


@define
class Iteration:

    dir_handle: DirHandle = field(converter=DirHandle)  # type: ignore
    iter_name: str = field(converter=str, kw_only=True)
    chainIDs: List[str] = field(kw_only=True, validator=validators.instance_of(list))
    resnames: List[str] = field(kw_only=True, validator=validators.instance_of(list))
    resSeqs: List[List[int]] = field(
        kw_only=True, validator=validators.instance_of(list)
    )  # Had to use validators.instance_of() instead of converter because mypy complains
    epoch_id: int = field(converter=int, init=False)
    complex: AbstractComplex = field(init=False)
    score_dir: DirHandle = field(converter=DirHandle, init=False)  # type: ignore
    scores: Dict[str, tuple] = field(init=False)
    mean_scores: Dict[str, float] = field(init=False)
    outside_box: bool = field(init=False, default=False)

    def __attrs_post_init__(self) -> None:
        try:
            self.epoch_id = int(Path(self.dir_handle).name.split("-")[0])
        except Exception as bad_iter_name:
            raise ValueError("Bad iteration name.") from bad_iter_name
        self.scores = {}
        self.mean_scores = {}

    def set_score(
        self, sf_name: str, scores: Sequence, log: Optional[logging.Logger] = None
    ) -> float:
        if not log:
            log = logging.getLogger("root")
        self.scores[sf_name] = tuple(scores)
        avg_score = mean(scores)
        self.mean_scores[sf_name] = avg_score

        std_score = stdev(scores)
        if abs(avg_score) < std_score:
            log.warning(
                f"{sf_name} score has a mean of {avg_score:.5f} and a std dev of {std_score:.5f}. "
                f"This is too much variance. You might want to check {self.epoch_id}-{self.iter_name}"
            )
        return self.mean_scores[sf_name]

    def write_down_scores(self):
        for sf_name, score in self.scores.items():
            with open(Path(self.score_dir, "scores_" + sf_name), "w") as f:
                for s in score:
                    f.write(str(round(s, 3)) + "\n")

    def read_scores(
        self, scoring_functions: Iterable, log: Optional[logging.Logger] = None
    ) -> None:
        if not log:
            log = logging.getLogger("root")
        if not getattr(self, "score_dir", None):
            try:
                self.score_dir = DirHandle(Path(self.dir_handle, "scoring"), make=False)
            except FileNotFoundError as e:
                raise FileNotFoundError(
                    f"read_scores(): {self.iter_name} has no scores."
                ) from e
        for SF in scoring_functions:
            try:
                scores_fn = Path(self.score_dir, "scores_" + SF)
                with open(scores_fn, "r") as f:
                    self.set_score(SF, [float(linea.strip()) for linea in f], log)
                if SF == "bluues":
                    scores_fn = Path(self.score_dir, "scores_" + "bmf")
                    with open(scores_fn, "r") as f:
                        self.set_score(
                            "bmf", [float(linea.strip()) for linea in f], log
                        )
            except FileNotFoundError as e:
                log.warning(f"{self.iter_name} was not scored with {SF}.")
                raise FileNotFoundError from e

    def __str__(self) -> str:
        return str(self.dir_handle)

    def __fspath__(self) -> str:
        return str(self.dir_handle)

    def __truediv__(self, key) -> Union[FileHandle, "DirHandle"]:
        return self.dir_handle.__truediv__(key)

    def __lt__(self, other) -> bool:
        return self.iter_name < other.iter_name


@define
class Epoch(MutableMapping):
    id: int = field(converter=int)
    iterations: Dict[str, Iteration] = field(kw_only=True)
    nvt_done: bool = field(converter=bool, default=False, kw_only=True)
    npt_done: bool = field(converter=bool, default=False, kw_only=True)
    top_iterations: Dict[str, Iteration] = field(init=False)
    mutated_positions: Set[int] = field(init=False)

    def set_top_iter(self, config_paths: Dict) -> None:
        # First, check if a "top_iterations" was set from a pickle tracking file.
        if "top_iterations" in config_paths:
            for top_it_full in config_paths["top_iterations"]:
                top_it_name = "-".join(Path(top_it_full).name.split("-")[1:])
                try:
                    self.top_iterations[top_it_name] = self.iterations[top_it_name]
                except Exception as e:
                    raise RuntimeError(
                        f"Could not set top iteration from pickle tracking file. Aborting."
                    ) from e
            return

        better_iters: PriorityQueue[Tuple[int, Iteration]] = PriorityQueue()

        for it in self.iterations.values():
            count = 0
            for other_it in self.iterations.values():
                if it == other_it:
                    continue
                # Allow the user to change scoring functions mid-run and only use
                # the subset present in both iterations.
                scoring_functions = set(it.scores.keys()) & set(other_it.scores.keys())
                count += sum(
                    [
                        other_it.mean_scores[SF] >= it.mean_scores[SF]
                        for SF in scoring_functions
                    ]
                )
            better_iters.put((-count, it))

        prev_count = 1
        while not better_iters.empty():
            # Remember, `count`, the priority, is negative.
            count, iter = better_iters.get()
            if count > prev_count:
                assert (
                    len(self.top_iterations) > 0
                ), f"Logical error. This can't happen."
                break
            prev_count = count
            self.top_iterations[iter.iter_name] = iter

    def backup(self) -> None:
        for it in self.iterations.values():
            orig = it.dir_handle.dir_path
            dest = Path(orig.parent, "bu_" + orig.name)
            if dest.exists():
                unique_str = "_".join(time.ctime().split()).replace(":", "_")
                dest = Path(orig.parent, f"bu_{unique_str}_{orig.name}")
            # TODO: won't be necessary to cast after 3.9 upgrade
            sh.move(str(orig), str(dest))

    def __getitem__(self, key) -> Iteration:
        return self.iterations[key]

    def __setitem__(self, key, value) -> None:
        self.iterations[key] = value

    def __delitem__(self, key) -> None:
        del self.iterations[key]

    def __iter__(self) -> Iterator:
        return self.iterations.__iter__()

    def __contains__(self, value: object) -> bool:
        return self.iterations.__contains__(value)

    def __len__(self) -> int:
        return len(self.iterations)


class WorkProject:
    """WorkProject main state holder of the run
    Args:
        config (Dict): dict with the yaml input configuration parsed by cli.py.
        start (bool): wether the working project is being run for the 1st time.
    Attributes:
            config (Dict): dict with the yaml input configuration parsed by cli.py.
            name (str): project's name
            dir_handle (DirHandle): working directory
            epochs (List[Epoch]): list of all the epochs generated on this run
            mdps (Dict[str, FileHandle]): dict with handles to the minimization,
                NVT and NPT mdps.
            scorers (Dict[str, AbstractScoringFunction]): dict with the callables to
                the scoring functions.
            last_mutated_position (Set[int]): last mutated position
            mutated_positions (Deque[Set[int]]): last N mutated positions
            mutated_aminoacids (Deque[Set[str]]): last N added amino acids. Unused for now.
            has_memory (bool = False): wether the protocol uses memory of mutations.
            has_failed_memory (bool = False): wether the protocol uses memory of failed mutations.
    Raises:
        e_start_cpx: Failure to build starting complex.
        e_restart_cpx: Failure to build a complex from input iterations.
        e_score: Failure to read scores from previous iterations.
        e_top_iter: Failure to set the top iterations from previous iterations.
        e: Other failures, like file issues or bad iteration names.
    """

    config: Dict
    name: str
    dir_handle: DirHandle
    epochs: List[Epoch]
    mdps: Dict[str, FileHandle]
    scorers: Dict[str, AbstractScoringFunction]
    mutated_positions: Deque[Set[int]] = deque(set())
    mutated_aminoacids: Deque[Set[str]] = deque(set())
    failed_mutated_positions: Deque[Set[int]] = deque(set())
    has_memory: bool = False
    has_failed_memory: bool = False

    def __init__(self, config: Dict, start: bool):
        self.config = config
        self.name = self.config["main"]["name"]
        self.epochs = []
        self.dir_handle = DirHandle(Path(self.config["paths"]["work"]), make=False)  # type: ignore
        log = logging.getLogger(self.name)

        if start:
            self.__start_work__(log)
        else:
            self.__restart_work__(log)

        self.get_mdps(self.config["paths"]["mdp"], self.config["md"]["mdp_names"])
        self.__add_scoring_functions__()
        self.__set_memory__()
        self.__set_failed_memory__()

        # TODO: brutto come la fame. checking `start` twice. Fix later.
        if start:
            log.info("Start of new Work Dir, not writing tracking info yet.")
        else:
            self.__track_project__(log)

    def __start_work__(self, log: logging.Logger):
        zero_epoch = Epoch(0, iterations={}, nvt_done=False, npt_done=False)
        for data_str in self.config["paths"]["input"]:

            # Check input PDB to create name and attributes for the starting iteration.
            input_path = Path(data_str)
            iter_name, chainIDs, resSeqs, resnames = self.__generate_iteration_ID__(
                input_path
            )

            iter_folder = "0-" + iter_name
            this_iter = Iteration(
                DirHandle(self.dir_handle.dir_path / iter_folder, make=True),
                iter_name=iter_name,
                chainIDs=chainIDs,
                resnames=resnames,
                resSeqs=resSeqs,
            )
            # Copy the input PDB into the iteration folder
            pdb_handle = FileHandle(input_path / (self.config["main"]["name"] + ".pdb"))
            copy_to(pdb_handle, this_iter.dir_handle)

            # set up complex
            try:
                this_iter.complex = GROComplex.from_pdb(
                    name=self.config["main"]["name"],
                    input_dir=this_iter.dir_handle,
                    target_chains=self.config["target"]["chainID"],
                    binder_chains=self.config["binder"]["chainID"],
                    md_config=self.config["md"],
                )
            except Exception as e_start_cpx:
                log.error(
                    f"Cannot generate starting complex from {pdb_handle}. Aborting."
                )
                raise e_start_cpx

            zero_epoch[iter_name] = this_iter

        zero_epoch.mutated_positions = set()
        # Finally, add the epoch 0.
        self.new_epoch(zero_epoch)

    def __restart_work__(self, log: logging.Logger):
        # Restart from input iterations, they should all have the same epoch number
        epoch_nbr = int(
            Path(self.config["paths"]["current_iterations"][0]).name.split("-")[0]
        )

        if "previous_iterations" in self.config["paths"]:
            prev_epoch = Epoch(
                epoch_nbr - 1, iterations={}, nvt_done=True, npt_done=True
            )
            prev_epoch.top_iterations = {}
            for iter_str in self.config["paths"]["previous_iterations"]:
                # Get `iter_name` from the input iteration dir
                # TODO: check this works when input iterations come from different epochs.
                iter_path = Path(iter_str)
                _, *name_by_chain = iter_path.name.split("-")
                iter_name = "-".join(name_by_chain)

                # `iter_name` comes from the input iteration.
                _, chainIDs, resSeqs, resnames = self.__generate_iteration_ID__(
                    iter_path
                )

                this_iter = Iteration(
                    DirHandle(iter_path, make=False),
                    iter_name=iter_name,
                    chainIDs=chainIDs,
                    resnames=resnames,
                    resSeqs=resSeqs,
                )
                try:
                    # Previous iterations should be fully ran.
                    this_iter.complex = GROComplex.from_complex(
                        name=self.config["main"]["prefix"]
                        + self.config["main"]["name"],
                        iter_path=iter_path,
                        target_chains=self.config["target"]["chainID"],
                        binder_chains=self.config["binder"]["chainID"],
                        ignore_cpt=False,
                    )
                except Exception as e_restart_cpx:
                    log.error(
                        f"Could not build complex from previous iteration: {iter_path}"
                    )
                    raise e_restart_cpx

                # Previous iterations should be fully scored.
                try:
                    this_iter.read_scores(self.config["scoring"]["functions"])
                except Exception as e_score:
                    log.error(f"Previous epoch not fully scored. Run in --score mode.")
                    raise e_score
                prev_epoch[iter_name] = this_iter

            try:
                prev_epoch.set_top_iter(self.config["paths"])
            except Exception as e_top_iter:
                log.error(f"Failed to set top iteration from: {prev_epoch.keys()}")
                raise e_top_iter

            top_itrs_str = " ; ".join(
                [
                    f"{iter.epoch_id}-{iter.iter_name}"
                    for iter in prev_epoch.top_iterations.values()
                ]
            )
            log.info(f"Previous epoch {epoch_nbr-1} top iterations: {top_itrs_str}")

            prev_epoch.mutated_positions = set()
            self.new_epoch(prev_epoch)

        current_epoch = Epoch(epoch_nbr, iterations={}, nvt_done=True, npt_done=True)
        for iter_str in self.config["paths"]["current_iterations"]:
            # Get `iter_name` from the input iteration dir
            iter_path = Path(iter_str)
            _, *name_by_chain = iter_path.name.split("-")
            iter_name = "-".join(name_by_chain)

            # `iter_name` comes from the input iteration.
            _, chainIDs, resSeqs, resnames = self.__generate_iteration_ID__(iter_path)

            this_iter = Iteration(
                DirHandle(iter_path, make=False),
                iter_name=iter_name,
                chainIDs=chainIDs,
                resnames=resnames,
                resSeqs=resSeqs,
            )
            # Create complex with coordinates and topology (should be zip)
            try:
                cpx_name = self.config["main"]["prefix"] + self.config["main"]["name"]
                this_iter.complex = GROComplex.from_complex(
                    name=cpx_name,
                    iter_path=iter_path,
                    target_chains=self.config["target"]["chainID"],
                    binder_chains=self.config["binder"]["chainID"],
                    ignore_cpt=False,
                )
            except Exception:
                try:
                    log.info(f"{iter_path.name} didn't finish its NPT MD.")
                    current_epoch.npt_done = False

                    cpx_name = "nvt_" + self.config["main"]["name"]
                    this_iter.complex = GROComplex.from_complex(
                        name=cpx_name,
                        iter_path=iter_path,
                        target_chains=self.config["target"]["chainID"],
                        binder_chains=self.config["binder"]["chainID"],
                        ignore_cpt=True,
                    )
                except:
                    try:
                        log.info(f"{iter_path.name} didn't finish its NVT MD.")
                        current_epoch.nvt_done = False
                        current_epoch.npt_done = False

                        cpx_name = self.config["main"]["name"]
                        this_iter.complex = GROComplex.from_complex(
                            name=cpx_name,
                            iter_path=iter_path,
                            target_chains=self.config["target"]["chainID"],
                            binder_chains=self.config["binder"]["chainID"],
                            ignore_cpt=True,
                        )
                    except Exception as e_restart_cpx:
                        log.error(
                            f"{iter_path} is in an invalid state. Cannot build complex from it."
                        )
                        raise e_restart_cpx

            current_epoch[iter_name] = this_iter

        current_epoch.mutated_positions = set(
            self.config["misc"]["epoch_mutated_positions"]
        )
        self.new_epoch(current_epoch)

    def __add_scoring_functions__(self) -> None:
        self.scorers = {}
        for sf in self.config["scoring"]["functions"]:
            self.scorers[str(sf)] = scoringfunctions[sf](
                self.config["paths"]["scoring_functions"],
                self.config["scoring"]["nprocs"],
            )

    def __get_mutating_resname__(
        self, pdb_path: Path, chainIDs: List, resSeqs: List
    ) -> List[str]:
        u = mda.Universe(str(pdb_path))

        resnames = []
        for chainID, list_resSeq in zip(chainIDs, resSeqs):
            ch_resnames = []
            ch_resids = []
            cadena = u.select_atoms(f"segid {chainID}")
            for resn, resi in {
                (atm.resname, atm.resnum)
                for atm in cadena.select_atoms(
                    " or ".join([f"resid {res}" for res in list_resSeq])
                )
            }:
                ch_resnames.append(resn)
                ch_resids.append(resi)
            resnames.append(
                "".join(
                    [
                        seq1(resn_3)
                        for resn_3 in np.array(ch_resnames)[np.argsort(ch_resids)]
                    ]
                )
            )
        return resnames

    def __generate_iteration_ID__(
        self, input_path: Path
    ) -> Tuple[str, List[str], List[List[int]], List[str]]:
        chainIDs = self.config["binder"]["mutating_chainID"]
        resSeqs = self.config["binder"]["mutating_resSeq"]
        try:
            pdb_path = input_path / (self.config["main"]["name"] + ".pdb")
            resnames = self.__get_mutating_resname__(pdb_path, chainIDs, resSeqs)
        except FileNotFoundError as e:
            # Give it another chance, we only need the PDB for topology info.
            pdb_path = input_path / (
                self.config["main"]["prefix"] + self.config["main"]["name"] + ".pdb"
            )
            resnames = self.__get_mutating_resname__(pdb_path, chainIDs, resSeqs)

        iter_name = "".join(
            [f"-{chainID}_{resname}" for chainID, resname in zip(chainIDs, resnames)]
        )

        # Drop the leading -:
        return iter_name[1:], chainIDs, resSeqs, resnames

    def __set_memory__(self):
        try:
            self.mutated_positions = deque(
                maxlen=self.config["protocol"]["memory_size"]
            )
            self.has_memory = True
        except KeyError:
            # User requested no memory
            self.mutated_positions = deque()
            return
        try:
            for set_of_positions in self.config["protocol"]["memory_positions"]:
                self.mutated_positions.appendleft(set(set_of_positions))
        except KeyError:
            self.mutated_positions.appendleft(set())
        # TODO: implement memory of amino acids?

    def __set_failed_memory__(self):
        try:
            self.failed_mutated_positions = deque(
                maxlen=self.config["protocol"]["failed_memory_size"]
            )
            self.has_failed_memory = True
        except KeyError:
            # User requested no memory
            self.failed_mutated_positions = deque()
            return
        try:
            for set_of_positions in self.config["protocol"]["failed_memory_positions"]:
                self.failed_mutated_positions.appendleft(set(set_of_positions))
        except KeyError:
            self.failed_mutated_positions.appendleft(set())

    def get_mem_aminoacids(self) -> Set[str]:
        set_of_aas: Set[str] = set()

        if not self.has_memory:
            return set_of_aas

        for mem in self.mutated_aminoacids:
            set_of_aas.update(mem)
        return set_of_aas

    def get_mem_positions(self) -> Set[int]:
        set_of_pos: Set[int] = set()

        if self.has_memory:
            for mem in self.mutated_positions:
                set_of_pos.update(mem)
        if self.has_failed_memory:
            for mem in self.failed_mutated_positions:
                set_of_pos.update(mem)
        return set_of_pos

    def get_mdps(self, dir_path: Path, mdp_names: Dict) -> None:
        self.mdps = {}
        for mdp_name, mdp_filename in mdp_names.items():
            mdp_path = Path(dir_path, mdp_filename)
            try:
                self.mdps[mdp_name] = FileHandle(mdp_path)
            except Exception as e:
                print(f"Could not get mdp file from: {mdp_path}.", flush=True)
                raise e

    def new_epoch(self, epoch: Epoch, log: Optional[logging.Logger] = None) -> None:
        self.epochs.append(epoch)
        self.__track_project__(log)

    def __track_project__(self, log: Optional[logging.Logger] = None) -> None:
        try:
            previous_iterations = [
                str(pre_it.dir_handle) for pre_it in self.epochs[-2].values()
            ]
            current_iterations = [
                str(cur_it.dir_handle) for cur_it in self.epochs[-1].values()
            ]
            top_iterations = [
                str(cur_it.dir_handle)
                for cur_it in self.epochs[-2].top_iterations.values()
            ]

            tracking = {
                "previous_iterations": previous_iterations,
                "current_iterations": current_iterations,
                "top_iterations": top_iterations,
                "epoch_mutated_positions": self.epochs[-1].mutated_positions,
                "memory_positions": self.mutated_positions,
                "failed_memory_positions": self.failed_mutated_positions,
            }
            assert (
                len(previous_iterations) > 0
                and len(current_iterations) > 0
                and len(top_iterations) > 0
            )

            # Back up tracking.pkl before writing.
            track = Path(self.dir_handle, "tracking.pkl")
            try:
                sh.move(track, Path(self.dir_handle, "bu_tracking.pkl"))
            except Exception:
                warn("No tracking.pkl file present. Writing initial one.")

            with open(track, "wb") as cur_file:
                pickle.dump(tracking, cur_file)

        except Exception as e:
            if log:
                log.warning(f"Not a full WorkProject yet, could not track it. {e}")

    def get_first_iter(self) -> Tuple[str, Iteration]:
        return next(iter(self.epochs[0].items()))


def set_logger(name: str, dir_path: Path) -> logging.Logger:
    logger = logging.getLogger(f"{name}")
    logger.setLevel(logging.INFO)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.WARNING)

    file_handler = logging.FileHandler(filename=f"{Path(dir_path, name)}.log")
    file_handler.setLevel(logging.INFO)

    formatter = logging.Formatter(
        "{asctime} - {levelname} - {message}",
        style="{",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    stream_handler.setFormatter(formatter)
    file_handler.setFormatter(formatter)

    # add Handlers to our logger
    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

    return logger
