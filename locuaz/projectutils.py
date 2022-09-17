from pathlib import Path
from attrs import define, field
from typing import Iterator, List, Sequence, Set, Dict, Tuple, Union, Deque
import collections
import logging
from statistics import mean, stdev
from collections import deque

from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

from fileutils import FileHandle, DirHandle, copy_to
from molecules import AbstractComplex, GROComplex
from abstractscoringfunction import AbstractScoringFunction
from scoringfunctions import scoringfunctions


@define
class Iteration:

    dir_handle: DirHandle = field(converter=DirHandle)  # type: ignore
    iter_name: str = field(converter=str, kw_only=True)
    chainIDs: List[str] = field(converter=list, kw_only=True)
    resnames: List[str] = field(converter=list, kw_only=True)
    resSeqs: List[List[int]] = field(converter=list, kw_only=True)
    complex: AbstractComplex = field(init=False)
    score_dir: DirHandle = field(converter=DirHandle, init=False)  # type: ignore
    scores: Dict[str, tuple] = field(init=False)
    mean_scores: Dict[str, float] = field(init=False)

    def __attrs_post_init__(self) -> None:
        self.scores = {}
        self.mean_scores = {}

    def set_score(self, sf_name: str, scores: Sequence) -> float:
        self.scores[sf_name] = tuple(scores)
        avg_score = mean(scores)
        self.mean_scores[sf_name] = avg_score

        std_score = stdev(scores)
        if abs(avg_score) < std_score:
            logging.warning(
                f"{sf_name} score has a mean of {avg_score} and a standard deviation "
                f"of {std_score}. This is too much variance. You might want to check this run. "
            )
        return self.mean_scores[sf_name]

    def write_down_scores(self):
        for sf_name, score in self.scores.items():
            with open(Path(self.score_dir, "scores_" + sf_name), "w") as f:
                for s in score:
                    f.write(str(round(s, 3)) + "\n")

    def read_scores(self, scoring_functions: Sequence) -> None:
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
                    self.set_score(SF, [float(linea.strip()) for linea in f])
            except FileNotFoundError as e:
                logging.warning(f"{self.iter_name} was not scored with {SF}.")

    def __str__(self) -> str:
        return str(self.dir_handle)

    def __fspath__(self) -> str:
        return str(self.dir_handle)

    def __truediv__(self, key) -> Union[FileHandle, "DirHandle"]:
        return self.dir_handle.__truediv__(key)

    def __lt__(self, other) -> bool:
        return self.iter_name < other.iter_name


@define
class Epoch(collections.abc.MutableMapping):
    id: int = field(converter=int)
    iterations: Dict[str, Iteration] = field()
    top_iterations: Dict[str, Iteration] = field(init=False)

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


@define
class WorkProject:
    """WorkProject:
    Args:
        config (Dict): dict with all the input configuration read from the input yaml.
    Attributes:
        root (Path): path to locuaz's root folder.
    """

    config: Dict = field()
    name: str = field(init=False)
    dir_handle: DirHandle = field(init=False)
    epochs: List[Epoch] = field(init=False)
    mdps: Dict[str, FileHandle] = field(init=False)
    history: Set[Tuple] = field(init=False)
    mutant_mtx: List[List[str]] = field(init=False)
    scorers: Dict[str, AbstractScoringFunction] = field(init=False)
    mutated_positions: Deque[Set[int]] = field(init=False)
    mutated_aminoacids: Deque[Set[str]] = field(init=False)
    has_memory: bool = field(init=False, default=False)

    def __attrs_post_init__(self):
        self.name = self.config["main"]["name"]
        self.epochs = []

        if "root" in self.config["paths"]:
            self.__start_work__()
        else:
            self.__restart_work__()

        self.get_mdps(self.config["paths"]["mdp"], self.config["md"]["mdp_names"])
        self.__add_scoring_functions__()
        self.__set_memory__()
        logging.basicConfig(
            filename=Path(self.dir_handle, f"{self.name}.log"), level=logging.INFO
        )

    def __start_work__(self):
        zero_epoch = Epoch(0, {})
        for data_str in self.config["paths"]["input"]:
            # First, create working dir
            self.dir_handle = DirHandle(
                Path(self.config["paths"]["root"], "work_dir"), make=True, force=True
            )

            # Check input PDB to create name and attributes for the starting iteration.
            input_path = Path(data_str)
            iter_name, chainIDs, resSeqs, resnames = self.__generate_iteration_ID__(
                input_path
            )

            self.history = set(iter_name)
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
            this_iter.complex = GROComplex.from_pdb(
                name=self.config["main"]["name"],
                input_dir=this_iter.dir_handle,
                target_chains=self.config["target"]["chainID"],
                binder_chains=self.config["binder"]["chainID"],
            )

            zero_epoch[iter_name] = this_iter

        # Finally, add the epoch 0.
        self.new_epoch(zero_epoch)

    def __restart_work__(self):
        # Restart from input iterations, but use the most recent epoch's number
        epoch_nbr = max(
            [
                int(Path(iteration_str).name.split("-")[0])
                for iteration_str in self.config["paths"]["current_iterations"]
            ]
        )
        if "previous_iterations" in self.config["paths"]:
            prev_epoch = Epoch(epoch_nbr - 1, {})
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
                # Create complex with coordinates (from npt run) and topology (should be zip)
                this_iter.complex = GROComplex.from_complex(
                    name=self.config["main"]["prefix"] + self.config["main"]["name"],
                    iter_path=iter_path,
                    target_chains=self.config["target"]["chainID"],
                    binder_chains=self.config["binder"]["chainID"],
                    ignore_cpt=False,
                )
                # Previous iterations should be complete.
                this_iter.read_scores(self.config["scoring"]["functions"])
                prev_epoch[iter_name] = this_iter
                prev_epoch.top_iterations[iter_name] = this_iter
            self.new_epoch(prev_epoch)

        current_epoch = Epoch(epoch_nbr, {})
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

            # Create complex with coordinates (from npt run) and topology (should be zip)
            this_iter.complex = GROComplex.from_complex(
                name=self.config["main"]["prefix"] + self.config["main"]["name"],
                iter_path=iter_path,
                target_chains=self.config["target"]["chainID"],
                binder_chains=self.config["binder"]["chainID"],
                ignore_cpt=False,
            )

            current_epoch[iter_name] = this_iter

        # Set up working dir
        self.dir_handle = DirHandle(iter_path.parent, make=False)
        self.new_epoch(current_epoch)

    # TODO: list the folders inside the working dir (assumed to be
    # `iter_path.parent`), discard those that aren't true iterationDirs
    # (in case the user had created some other folder), get the iternames
    # of those folder and add them to the history set.
    # self.history = all iternames in iter_path.parent
    def __add_scoring_functions__(self) -> None:
        first_iter = tuple(self.epochs[-1].iterations.values())[0]
        target_chains = tuple(first_iter.complex.top.target_chains)
        binder_chains = tuple(first_iter.complex.top.binder_chains)

        self.scorers = {}
        for sf in self.config["scoring"]["functions"]:
            self.scorers[str(sf)] = scoringfunctions[sf](
                self.config["paths"]["scoring_functions"],
                self.config["scoring"]["nprocs"],
            )

    def __get_mutating_resname__(
        self, pdb_path: Path, chainIDs: List, resSeqs: List
    ) -> List[str]:
        parsero = PDBParser(QUIET=True)
        pdb = parsero.get_structure("UNK", pdb_path)

        resnames = []
        for (chainID, resSeq) in zip(chainIDs, resSeqs):
            resname_by_chain = ""
            for chain in pdb.get_chains():
                if chain.get_id() == chainID:
                    for residue in chain.get_residues():
                        if residue.get_id()[1] in resSeq:
                            resname_by_chain += seq1(residue.get_resname())
            if len(resname_by_chain) == 0:
                raise ValueError(
                    f"Could not get resnames for the mutating resSeq: {resSeq} "
                    f"and chainID: {chainID}. mutating_resSeq may be wrong. Aborting."
                )
            resnames.append(resname_by_chain)

        self.config["binder"]["mutating_resnames"] = resnames
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
            return
        try:
            for set_of_positions in self.config["protocol"]["memory_positions"]:
                self.mem_positions(set_of_positions)
        except KeyError:
            self.mem_positions([])
        try:
            self.mutated_aminoacids = deque(
                maxlen=self.config["protocol"]["memory_size"]
            )
            for set_of_aas in self.config["protocol"]["memory_aminoacids"]:
                self.mem_aminoacids(set_of_aas)
            self.has_memory = True
        except KeyError:
            self.mem_aminoacids([])

    def mem_aminoacids(self, aa_set: Sequence) -> None:
        self.mutated_aminoacids.append(set(aa_set))

    def mem_positions(self, pos_set: Sequence) -> None:
        self.mutated_positions.append(set(pos_set))

    def get_mem_aminoacids(self) -> Set[str]:
        set_of_aas: Set[str] = set()

        if not self.has_memory:
            return set_of_aas

        for mem in self.mutated_aminoacids:
            set_of_aas.update(mem)
        return set_of_aas

    def get_mem_positions(self) -> Set[int]:
        set_of_pos: Set[int] = set()

        if not self.has_memory:
            return set_of_pos

        for mem in self.mutated_positions:
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

    def new_epoch(self, epoch: Epoch) -> None:
        self.epochs.append(epoch)

    def get_first_iter(self) -> Tuple[str, Iteration]:
        return next(iter(self.epochs[0].items()))
