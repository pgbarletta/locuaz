import logging
import os
import pickle
import shutil as sh
import time
from collections import deque
from collections.abc import MutableMapping
# TODO: replace own pairwise with itertools' on 3.10
# from itertools import pairwise
from itertools import tee, zip_longest
from pathlib import Path
from queue import PriorityQueue
from typing import (
    Iterable,
    Iterator,
    List,
    Set,
    Dict,
    Tuple,
    Union,
    Deque,
    Optional,
)
from warnings import warn
from numpy.typing import NDArray

import MDAnalysis as mda
import numpy as np
from Bio.SeqUtils import seq1
from attrs import define, field, validators
import matplotlib.pyplot as plt
import networkx as nx

from locuaz.abstractscoringfunction import AbstractScoringFunction
from locuaz.complex import AbstractComplex, GROComplex
from locuaz.fileutils import FileHandle, DirHandle, copy_to
from locuaz.primitives import AA_MAP, my_seq1
from locuaz.scoringfunctions import scoringfunctions
from locuaz.interface import get_freesasa_residues
from locuaz.mutation import Mutation


# TODO: replace own pairwise with itertools' on 3.10
def pairwise(iterable):
    # pairwise('ABCDEFG') --> AB BC CD DE EF FG
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


@define
class Iteration:
    dir_handle: DirHandle = field(converter=DirHandle)  # type: ignore
    iter_name: str = field(converter=str, kw_only=True)
    chainIDs: List[str] = field(kw_only=True, validator=validators.instance_of(list))
    resnames: List[str] = field(kw_only=True, validator=validators.instance_of(list))
    resSeqs: List[List[int]] = field(
        kw_only=True, validator=validators.instance_of(list)
    )  # Had to use validators.instance_of() instead of converter because mypy complains
    parent: Optional["Iteration"] = field(kw_only=True, default=None)
    mutation: Optional[Mutation] = field(kw_only=True, default=None)
    epoch_id: int = field(converter=int, init=False)
    complex: Union[AbstractComplex, GROComplex] = field(init=False)
    score_dir: DirHandle = field(converter=DirHandle, init=False)  # type: ignore
    scores: Dict[str, tuple] = field(init=False)
    mean_scores: Dict[str, float] = field(init=False)
    stats: Dict[str, NDArray] = field(init=False)
    outside_box: bool = field(init=False, default=False)

    def __attrs_post_init__(self) -> None:
        try:
            self.epoch_id = int(Path(self.dir_handle).name.split("-")[0])
        except Exception as bad_iter_name:
            raise ValueError("Bad iteration name.") from bad_iter_name
        self.scores = {}
        self.mean_scores = {}
        self.stats = {}

    def set_score(
            self, sf_name: str, scores: Iterable, log: Optional[logging.Logger] = None
    ) -> float:
        if not log:
            log = logging.getLogger("root")
        self.scores[sf_name] = tuple(scores)
        avg_score: float = np.mean(scores)
        self.mean_scores[sf_name] = avg_score

        std_score = np.std(scores)
        if abs(avg_score) < std_score:
            log.warning(
                f"{sf_name} score has a mean of {avg_score:.5f} and a std dev of {std_score:.5f}. "
                f"This is too much variance. You might want to check {self.epoch_id}-{self.iter_name}"
            )
        return self.mean_scores[sf_name]

    def set_stat(
            self, stat_name: str, stat: NDArray[float], log: Optional[logging.Logger] = None
    ) -> float:
        self.stats[stat_name] = stat
        avg_stat: float = np.mean(stat)
        if log:
            log.info(f"'{stat_name}': {avg_stat:.3f} +/- {np.std(stat):.3f}")
        return avg_stat

    def write_down_scores(self):
        for sf_name, score in self.scores.items():
            with open(Path(self.score_dir, "scores_" + sf_name), "w") as f:
                for s in score:
                    f.write(str(round(s, 3)) + "\n")

    def write_down_stats(self):
        for stat_name, result in self.stats.items():
            with open(Path(self.score_dir, f"stat_{stat_name}"), "w") as f:
                for val in result:
                    f.write(str(round(val, 3)) + "\n")

    def read_scores(
            self, scoring_functions: Iterable, log: Optional[logging.Logger] = None
    ) -> bool:
        if not log:
            log = logging.getLogger("root")
        if not getattr(self, "score_dir", None):
            try:
                self.score_dir = DirHandle(Path(self.dir_handle, "scoring"), make=False)
            except FileNotFoundError as e:
                # No scoring folder
                return False
                # raise FileNotFoundError(
                #     f"read_scores(): {self.iter_name} has no scores."
                # ) from e
        for SF in scoring_functions:
            try:
                scores_fn = Path(self.score_dir, "scores_" + SF)
                with open(scores_fn, "r") as f:
                    self.set_score(SF, [float(linea.strip()) for linea in f], log)
            except FileNotFoundError as e:
                log.warning(f"{self.epoch_id}-{self.iter_name} was not scored with {SF}.")
                return False
        return True

    def mutation_str(self) -> Optional[str]:
        if self.mutation:
            return self.mutation.to_str()
        else:
            return None

    def full_name(self) -> str:
        return f"{self.epoch_id}-{self.iter_name}"

    def __str__(self) -> str:
        return str(self.dir_handle)

    def __fspath__(self) -> str:
        return str(self.dir_handle)

    def __truediv__(self, key) -> Union[FileHandle, "DirHandle"]:
        return self.dir_handle.__truediv__(key)

    # noinspection PyDataclass
    def __lt__(self, other: "Iteration") -> bool:
        """
        Iterations will be added to a heapq during pruning. In case they have the same score,
        iterations will be compared against each other, so a comparison function is needed, even
        though we don't care about the order of iterations with the same score.
        Args:
            other (Iteration): another iteration object.
        Returns:
            bool: first iteration in alphabetical order according to its name.
        """
        return self.iter_name < other.iter_name

    def generate_name_resname(self, mutation: Mutation) -> Tuple[str, List[str]]:
        """
        returns the name of the iteration that would be generated by the input mutation.
        It also returns a list with the resnames (1-letter names) of each mutated sequence (CDR),
        also including the one that would be mutated.
        Parameters
        ----------
        mutation : Mutation

        Returns
        -------
        iteration_name and resnames: Tuple[str, List[str]]
        """
        iter_name = ""
        new_iteration_resnames = []
        for idx, (chainID, resname) in enumerate(zip(self.chainIDs, self.resnames)):
            if idx == mutation.chainID_idx:
                # This is the mutated chainID
                new_resname = resname[: mutation.resSeq_idx] + mutation.new_aa + resname[mutation.resSeq_idx + 1:]
            else:
                # This one remains the same
                new_resname = "".join([residue for residue in resname])
            iter_name += f"-{chainID}_{new_resname}"
            new_iteration_resnames.append(new_resname)

        # Drop the leading '-'
        return iter_name[1:], new_iteration_resnames


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
            count, iteration = better_iters.get()
            if count > prev_count:
                assert (
                        len(self.top_iterations) > 0
                ), f"Logical error. This can't happen."
                break
            prev_count = count
            self.top_iterations[iteration.iter_name] = iteration

    def backup(self) -> None:
        for it in self.iterations.values():
            orig = it.dir_handle.dir_path
            dest = Path(orig.parent, "bu_" + orig.name)
            if dest.exists():
                unique_str = "_".join(time.ctime().split()).replace(":", "_")
                dest = Path(orig.parent, f"bu_{unique_str}_{orig.name}")
            # TODO: won't be necessary to cast after 3.9 upgrade
            sh.move(orig, dest)

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
    config: Dict
    name: str
    dir_handle: DirHandle
    epochs: List[Epoch]
    mdps: Dict[str, FileHandle]
    scorers: Dict[str, AbstractScoringFunction] = {}
    project_dag: nx.DiGraph = nx.DiGraph()
    project_mut_dag: nx.DiGraph = nx.DiGraph()
    mutated_positions: Deque[Set[int]] = deque(set())
    mutated_aminoacids: Deque[Set[str]] = deque(set())
    failed_mutated_positions: Deque[Set[int]] = deque(set())
    has_memory: bool = False
    has_failed_memory: bool = False
    tleap_dir: Optional[DirHandle] = None

    def __init__(self, config: Dict, start: bool):
        self.config = config
        self.name = self.config["main"]["name"]
        self.epochs = []
        self.dir_handle = DirHandle(Path(self.config["paths"]["work"]), make=False)  # type: ignore
        log = logging.getLogger(self.name)
        self.__init_tleap__()

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
        starting_epoch = self.config["main"]["starting_epoch"]
        zero_epoch = Epoch(starting_epoch, iterations={}, nvt_done=False, npt_done=False)
        for data_str in self.config["paths"]["input"]:

            # Check input PDB to create name and attributes for the starting iteration.
            input_path = Path(data_str)
            try:
                self.__check_input_pdb__(input_path)
            except (Exception,) as e:
                # Remove the working dir, so it can restart easily again.
                sh.rmtree(Path(self.config["paths"]["work"]))
                raise e
            iter_name, chainIDs, resSeqs, resnames = self.__generate_iteration_ID__(input_path)

            iter_path = Path(self.dir_handle, f"{starting_epoch}-{iter_name}")
            assert not (iter_path.is_dir() or iter_path.is_file()), f"{iter_path} exists. Wrong input path."
            this_iter = Iteration(
                DirHandle(iter_path, make=True),
                iter_name=iter_name,
                chainIDs=chainIDs,
                resnames=resnames,
                resSeqs=resSeqs,
            )
            # Copy the input PDB into the iteration folder
            pdb_handle = FileHandle(input_path / (self.config["main"]["name"] + ".pdb"))
            copy_to(pdb_handle, Path(this_iter.dir_handle))

            # Copy tleap files, if necessary
            self.get_tleap_into_iter(Path(this_iter.dir_handle))

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
        # Finally, add the initial epoch.
        self.new_epoch(zero_epoch)

    def __restart_work__(self, log: logging.Logger):
        # Restart from input iterations, they should all have the same epoch number
        epoch_nbr = int(Path(self.config["paths"]["current_iterations"][0]).name.split("-")[0])

        if "previous_iterations" in self.config["paths"]:
            prev_epoch = Epoch(epoch_nbr - 1, iterations={}, nvt_done=True, npt_done=True)
            prev_epoch.top_iterations = {}
            for iter_str in self.config["paths"]["previous_iterations"]:
                # Get `iter_name` from the input iteration dir
                iter_name, iter_path, this_iter = self.iteration_from_str(iter_str)
                try:
                    # Previous iterations should be fully ran.
                    this_iter.complex = GROComplex.from_complex(
                        name=self.config["main"]["prefix"] + self.config["main"]["name"],
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

                # Finally, add its mutation, if present.
                try:
                    this_iter.mutation = self.config["mutations"][iter_name]
                except KeyError:
                    this_iter.mutation = None
                # Store it in the epoch.
                prev_epoch[iter_name] = this_iter

            try:
                prev_epoch.set_top_iter(self.config["paths"])
            except Exception as e_top_iter:
                log.error(f"Failed to set top iteration from: {prev_epoch.keys()}")
                raise e_top_iter

            top_itrs_str = " ; ".join(
                [
                    f"{iteration.epoch_id}-{iteration.iter_name}"
                    for iteration in prev_epoch.top_iterations.values()
                ]
            )
            log.info(f"Previous epoch {epoch_nbr - 1} top iterations: {top_itrs_str}")

            prev_epoch.mutated_positions = set()
            self.new_epoch(prev_epoch)

        current_epoch = Epoch(epoch_nbr, iterations={}, nvt_done=True, npt_done=True)
        for iter_str in self.config["paths"]["current_iterations"]:
            # Get `iter_name` from the input iteration dir
            iter_name, iter_path, this_iter = self.iteration_from_str(iter_str)
            # Create complex with coordinates and topology (should be .zip)
            try:
                cpx_name = self.config["main"]["prefix"] + self.config["main"]["name"]
                this_iter.complex = GROComplex.from_complex(
                    name=cpx_name,
                    iter_path=iter_path,
                    target_chains=self.config["target"]["chainID"],
                    binder_chains=self.config["binder"]["chainID"],
                    ignore_cpt=False,
                )
            except RuntimeError:
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
                except RuntimeError:
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
            # Finally, add its mutation, if present.
            try:
                this_iter.mutation = self.config["mutations"][iter_name]
            except KeyError:
                this_iter.mutation = None
            # Store it in the epoch.
            current_epoch[iter_name] = this_iter

        current_epoch.mutated_positions = set(self.config["misc"]["epoch_mutated_positions"])

        self.new_epoch(current_epoch)
        self.project_dag = self.config["project_dag"]
        self.project_mut_dag = self.config["project_mut_dag"]

    def iteration_from_str(self, iter_str: str) -> Tuple[str, Path, Iteration]:
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
        return iter_name, iter_path, this_iter

    def __add_scoring_functions__(self) -> None:
        for sf in self.config["scoring"]["functions"]:
            self.scorers[str(sf)] = scoringfunctions[sf](
                self.config["paths"]["scoring_functions"],
                nthreads=self.config["scoring"]["nthreads"],
                mpiprocs=self.config["scoring"]["mpiprocs"],
            )

    @staticmethod
    def __get_mutating_resname__(
            pdb_path: Path, chainIDs: List, resSeqs: List
    ) -> List[str]:

        u = mda.Universe(str(pdb_path))
        resnames = []
        for chainID, list_resSeq in zip(chainIDs, resSeqs):
            ch_resnames = []
            ch_resids = []
            cadena = u.select_atoms(f"segid {chainID}")
            for resn, resi in {(atm.resname, atm.resnum)
                               for atm in cadena.select_atoms(" or ".join([f"resid {res}" for res in list_resSeq]))}:
                ch_resnames.append(my_seq1(resn))
                ch_resids.append(resi)

            resnames.append(
                "".join([resn_1 for resn_1 in np.array(ch_resnames)[np.argsort(ch_resids)]])
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
        except FileNotFoundError:
            # Give it another chance, we only need the PDB for topology info.
            pdb_path = input_path / (
                    self.config["main"]["prefix"] + self.config["main"]["name"] + ".pdb"
            )
            resnames = self.__get_mutating_resname__(pdb_path, chainIDs, resSeqs)

        iter_name = "-".join(
            [f"{chainID}_{resname}" for chainID, resname in zip(chainIDs, resnames)]
        )
        return iter_name, chainIDs, resSeqs, resnames

    def __check_input_pdb__(self, input_path: Path) -> None:

        # Check the chainIDs:
        pdb_path = input_path / (self.config["main"]["name"] + ".pdb")
        u = mda.Universe(str(pdb_path))
        segids = [s.segid for s in u.segments]  # type: ignore
        assert len(segids) > 2, "Too few segments in the input PDB. " \
                                "There should be at least 3 (target+binder+solvent)."

        chainIDs = self.config["target"]["chainID"] + self.config["binder"]["chainID"]
        for segid, chainID in zip_longest(segids, chainIDs):
            if chainID:
                assert (segid == chainID), f"PDBs chainIDs ({segids}) and input target-binder chainIDs " \
                                           f"({chainIDs}) from the config should have the same ordering and the target chains should go first."
            else:
                assert (segid == '' or segid == 'X'), f"There're unaccounted segments. {segid} has to either be, " \
                                                      f"target or binder. If it's solvent or ions, then its segid " \
                                                      "and chainID should be empty ('') or 'X'."

        # Check amino acids
        all_residues = "protein or " + " or ".join(
            [f"segid {chainID}" for chainID in (self.config["target"]["chainID"] + self.config["binder"]["chainID"])])
        prot = u.atoms.select_atoms(all_residues)  # type: ignore
        aas = {res.resname for res in prot.residues}
        all_aas = set(AA_MAP.keys())
        nonstandard_residues = aas - all_aas
        if len(nonstandard_residues) != 0:
            warn(f"Unrecognized residues: {nonstandard_residues}. Make sure you can build "
                 "a topology for them and that they are not necessary for scoring.")

            allowed_nonstandard_residues = set(self.config["scoring"]["allowed_nonstandard_residues"])
            if not nonstandard_residues.issubset(allowed_nonstandard_residues):
                warn(f"Watch out, {nonstandard_residues - allowed_nonstandard_residues} were not included as "
                     "'allowed_nonstandard_residues', make sure you don't need them for scoring.")

        # Check solvent
        solvent_sel = "resname WAT" if self.config["md"]["use_tleap"] else "resname SOL"
        wat = u.select_atoms(solvent_sel)
        assert len(wat), f"Input PDB lacks solvent: {solvent_sel}."

        # Make sure mutating residues are present.
        mutating_chainID = self.config["binder"]["mutating_chainID"]
        mutating_resSeq = self.config["binder"]["mutating_resSeq"]
        mutating_resname = self.config["binder"]["mutating_resname"]

        # And also check against freesasa
        freesasa_resis = get_freesasa_residues(pdb_path, mutating_chainID)

        resis_not_present = set()
        resis_not_present_freesasa = set()
        correct_resis_resSeq = []
        correct_resis_resname = []
        for chainID, resSeqs, resnames in zip(mutating_chainID, mutating_resSeq, mutating_resname):
            for chain in u.segments:  # type: ignore
                if chain.segid == chainID:
                    selected_resis = {(resSeq, resname) for resSeq, resname in zip(resSeqs, resnames)}
                    available_resis = {(residue.resnum, seq1(residue.resname)) for residue in chain.residues}

                    bad_resis = selected_resis - available_resis
                    resis_not_present.update(bad_resis)
                    if len(bad_resis) > 0:
                        correct_resis = [(residue.resnum, seq1(residue.resname)) for residue in chain.residues if
                                         residue.resnum in resSeqs]
                        correct_resis_resSeq += [res[0] for res in correct_resis]
                        correct_resis_resname += [res[1] for res in correct_resis]

                    resis_not_present_freesasa.update(selected_resis - freesasa_resis)

        if len(resis_not_present) > 0:
            raise ValueError(
                f"The following residues are not present in the input PDB: {resis_not_present}.\n"
                f"This is the 'mutating_resname' for the input 'mutating_resSeq' "
                f"'{correct_resis_resSeq}': {correct_resis_resname}\n"
                "Don't just copy them, make sure your trying to mutate the right residues.")

        if len(resis_not_present_freesasa) > 0:
            raise ValueError(
                f"BUG: the following residues are not present according to freesasa: {resis_not_present_freesasa}.\n"
                f"These are the (resSeq, resname) pairs for chains {mutating_chainID} "
                f"according to freesasa: {sorted(freesasa_resis, key=lambda x: x[0])}")

        # Check residue numbering
        if self.config["md"]["use_tleap"]:
            # Check if resSeq is continuous.
            for i, j in pairwise(range(u.segments.n_segments)):  # type: ignore
                pre_seg = u.segments[i]  # type: ignore
                seg = u.segments[j]  # type: ignore
                if pre_seg.segid in chainIDs and seg.segid in chainIDs:
                    pre_resSeq = pre_seg.residues[-1].resnum
                    resSeq = seg.residues[0].resnum
                    assert (
                                   resSeq - pre_resSeq
                           ) == 1, f"Non-continuous resSeq ({pre_resSeq}, {resSeq}) between {pre_seg} and {seg}. "
                    "It should be continuous, according to Amber specifications. "
                    "'mutating_resSeq' should follow this same convention."
        else:
            # Check if each protein chain begins at resnum 1
            for segment in u.segments:  # type: ignore
                if segment.segid in chainIDs:
                    first_resSeq = segment.residues[0].resnum
                    assert (
                            first_resSeq == 1
                    ), f"First resSeq from chain {segment.segid} is {first_resSeq}. "
                    "It should be 1, according to GROMACS specifications. "
                    "'mutating_resSeq' should follow this same convention."

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

    def __init_tleap__(self) -> None:
        if self.config["md"]["use_tleap"]:
            self.tleap_dir = DirHandle(Path(self.config["paths"]["tleap"]))

    def get_tleap_into_iter(self, iter_dir: Path) -> None:
        if self.config["md"]["use_tleap"]:
            for file in os.listdir(Path(self.tleap_dir)):  # type: ignore
                sh.copy(Path(self.tleap_dir, file), Path(iter_dir))  # type: ignore

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
        self.update_dags()
        self.__track_project__(log)

    def update_dags(self) -> None:

        for iteration in self.epochs[-1].values():
            # Iteration DAG:
            try:
                self.project_dag.add_edge(iteration.parent.full_name(), iteration.full_name())
            except AttributeError:
                # This iteration is on the 0 epoch, it has no parent.
                self.project_dag.add_node(iteration.full_name())

            # Mutational DAG:
            try:
                self.project_mut_dag.add_edge(iteration.parent.mutation_str(), iteration.mutation_str())
            except ValueError:
                # The parent iteration is on the 0 epoch and has no mutation.
                self.project_mut_dag.add_edge("None", iteration.mutation_str())
            except AttributeError:
                # This iteration is on the 0 epoch, it has no parent and no mutation.
                self.project_mut_dag.add_node("None")

    def __track_project__(self, log: Optional[logging.Logger] = None) -> None:
        try:
            previous_iterations: List[str] = []
            mutations: Dict[str, str] = {}
            for iter_name, iteration in self.epochs[-2].items():
                previous_iterations.append(str(iteration.dir_handle))
                mutations[iter_name] = iteration.mutation

            current_iterations: List[str] = []
            for iter_name, iteration in self.epochs[-1].items():
                current_iterations.append(str(iteration.dir_handle))
                mutations[iter_name] = iteration.mutation

            top_iterations: List[str] = []
            for iter_name, iteration in self.epochs[-2].top_iterations.items():
                top_iterations.append(str(iteration.dir_handle))
                mutations[iter_name] = iteration.mutation

            tracking = {
                "previous_iterations": previous_iterations,
                "current_iterations": current_iterations,
                "top_iterations": top_iterations,
                "mutations": mutations,
                "epoch_mutated_positions": self.epochs[-1].mutated_positions,
                "memory_positions": self.mutated_positions,
                "failed_memory_positions": self.failed_mutated_positions,
                "project_dag": self.project_dag,
                "project_mut_dag": self.project_mut_dag,
            }
            assert (len(previous_iterations) > 0 and len(current_iterations) > 0 and len(top_iterations) > 0)

            # Back up tracking.pkl before writing.
            track = Path(self.dir_handle, "tracking.pkl")
            try:
                unique_str = "_".join(time.ctime().split()).replace(":", "_")
                sh.move(track, Path(self.dir_handle, f"{unique_str}_tracking.pkl"))
            except (FileNotFoundError, OSError):
                pass

            with open(track, "wb") as cur_file:
                pickle.dump(tracking, cur_file)

            self.__draw_dags__(Path(self.dir_handle, "dags.png"))
        except Exception as e:
            if log:
                log.warning(f"Not a full WorkProject yet, could not track it.")

    def get_first_iter(self) -> Tuple[str, Iteration]:
        return next(iter(self.epochs[0].items()))

    def __draw_dags__(self, out_path: Path) -> None:
        # Try to set a nice plot size
        w = max(dict(self.project_dag.degree()).values()) * 2 + 4
        h = (nx.dag_longest_path_length(self.project_dag) + 1) * 4 + 2
        _, axes = plt.subplots(1, 2, figsize=(w, h))
        # _, axes = plt.subplots(1, 2)

        ##### Iterations graph
        plt.sca(axes[0])
        # For a better display of the iterations names.
        label_remap = {node: node.replace('-', '\n') for node in self.project_dag.nodes}
        dag = nx.relabel_nodes(self.project_dag, label_remap)
        # Adjuste node and plot sizes, so it fits nicely, hopefully.
        one_label = next(iter(label_remap.keys()))
        node_size = max([len(piece) for piece in one_label.split('\n')]) * 400
        # Plot
        pos = nx.drawing.nx_agraph.graphviz_layout(dag, prog="dot")
        nx.draw(dag, pos, with_labels=True, node_size=node_size, node_color="lightblue", font_size=10)

        ##### Mutations graph
        plt.sca(axes[1])
        # Plot
        pos = nx.drawing.nx_agraph.graphviz_layout(self.project_mut_dag, prog="dot")
        nx.draw(self.project_mut_dag, pos, with_labels=True, node_size=4000, node_color="pink", font_size=10)

        # Remove axes and save.
        axes[0].set_axis_off()
        axes[1].set_axis_off()
        plt.savefig(out_path)

    @staticmethod
    def draw_dag(dag: nx.Graph, out_path: Path, *, reformat: bool = False) -> None:
        node_size = 4000
        w = max(dict(dag.degree()).values()) * 2 + 2
        h = (nx.dag_longest_path_length(dag) + 1) * 4 + 2
        plt.figure(figsize=(w, h))

        if reformat:
            # For a better display of the iterations names.
            label_remap = {node: node.replace('-', '\n') for node in dag.nodes}
            dag = nx.relabel_nodes(dag, label_remap)

            # Adjuste node and plot sizes, so it fits nicely, hopefully.
            one_label = next(iter(label_remap.keys()))
            node_size = max([len(piece) for piece in one_label.split('\n')]) * 500

        # Plot and save
        pos = nx.drawing.nx_agraph.graphviz_layout(dag, prog="dot")
        nx.draw(dag, pos, with_labels=True, node_size=node_size, node_color="lightblue", font_size=10)
        plt.savefig(out_path)


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
