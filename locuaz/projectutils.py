from pathlib import Path
import shutil as sh
from attrs import define, field
from typing import Iterator, List, Sequence, Set, Dict, Tuple
import collections
from abc import ABCMeta
import logging
from statistics import mean
from fileutils import FileHandle, DirHandle, copy_to
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

from primitives import launch_biobb
from molecules import AbstractComplex, GROComplex, PDBStructure
from abstractscoringfunction import AbstractScoringFunction
from scoringfunctions import *


@define
class AbstractDir(metaclass=ABCMeta):
    dir_handle: DirHandle = field(converter=DirHandle)  # type: ignore
    name: str = field(init=False)

    def __attrs_post_init__(self):
        self.name = self.dir_handle.dir_path.name


@define
class Iteration(AbstractDir):

    iter_name: str = field(converter=str, kw_only=True)
    chainIDs: List[str] = field(converter=list, kw_only=True)
    resnames: List[str] = field(converter=list, kw_only=True)
    resSeqs: List[List[int]] = field(converter=list, kw_only=True)
    complex: AbstractComplex = field(init=False)
    score_dir: DirHandle = field(converter=DirHandle, init=False)  # type: ignore
    scores: Dict[str, tuple] = field(init=False)
    mean_scores: Dict[str, float] = field(init=False)
    # TODO: this won't work if the complex has nucleic acids or glyco stuff.
    nonwat_pdb: PDBStructure = field(init=False)
    wation_pdb: PDBStructure = field(init=False)

    def __attrs_post_init__(self) -> None:
        self.scores = {}
        self.mean_scores = {}

    def set_score(self, sf_name: str, scores: Sequence) -> float:
        self.scores[sf_name] = tuple(scores)
        self.mean_scores[sf_name] = mean(scores)
        return self.mean_scores[sf_name]

    def write_down_scores(self):
        print(f"--- PARCE info: writing output scores.")

        for sf_name, score in self.scores.items():
            with open(Path(self.score_dir, "scores_" + sf_name), "w") as f:
                for s in score:
                    f.write(str(round(s, 3)) + "\n")


@define
class Epoch(collections.abc.MutableMapping):
    id: int = field(converter=int)
    iterations: Dict[str, Iteration] = field()
    top_iterations: Tuple[str] = field(init=False)

    def __getitem__(self, key):
        return self.iterations[key]

    def __setitem__(self, key, value):
        self.iterations[key] = value

    def __delitem__(self, key):
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

    def __attrs_post_init__(self):
        self.name = self.config["main"]["name"]
        self.epochs = []

        if self.config["main"]["mode"] == "start":
            self.__start_work__()
        elif self.config["main"]["mode"] in ("restart", "score"):
            self.__restart_work__()
        else:
            raise RuntimeError("Horrible bug, this shouldn't happen.")

        # MDP files
        self.update_all_mdps(
            self.config["paths"]["data"], self.config["md"]["mdp_names"]
        )
        self.__add_scoring_functions__()

    def __start_work__(self):
        # First, create working dir
        self.dir_handle = DirHandle(
            Path(self.config["paths"]["root"], "work_dir"), make=True, force=True
        )
        # Set up the logger:
        logging.basicConfig(
            filename=Path(self.dir_handle, "root_log"), level=logging.INFO
        )
        logging.info(f"Starting from work project: {self.dir_handle}")

        # Then, create the dir for the 0 iteration (original input files) and set it up
        input_path = Path(self.config["paths"]["input"][0])
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

        # Finally, add the epoch 0.
        self.new_epoch(Epoch(0, {iter_name: this_iter}))

    def __restart_work__(self):
        current_epoch = Epoch(self.config["epoch_nbr"], {})

        for iter_str in self.config["paths"]["input"]:
            # Get `iter_name` from the input iteration dir
            iter_path = Path(iter_str)
            _, *name_by_chain = iter_path.name.split("-")
            iter_name = "-".join(name_by_chain)

            # Set up working dir
            self.dir_handle = DirHandle(iter_path.parent, make=False)

            # Set up the logger:
            logging.basicConfig(
                filename=Path(self.dir_handle, "root_log"), level=logging.INFO
            )
            logging.info(f"Restarting from work project: {self.dir_handle}")

            # `iter_name` is already defined.
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
            )

            current_epoch[iter_name] = this_iter

        # Finally, add the corresponding epoch.
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
            sf_class = eval(sf)

            self.scorers[str(sf) + "_scorer"] = sf_class(
                self.config["paths"]["scoring_functions"],
                self.config["scoring"]["nprocs"],
                target_chains=target_chains,
                binder_chains=binder_chains,
            )

    def __get_mutating_resname__(
        self, pdb_path: Path, chainIDs: List, resSeqs: List
    ) -> List[str]:
        parsero = PDBParser(QUIET=True)
        pdb = parsero.get_structure("1ppg", pdb_path)

        resnames = []
        for (chainID, resSeq) in zip(chainIDs, resSeqs):
            resname_by_chain = ""
            for chain in pdb.get_chains():
                if chain.get_id() == chainID:
                    for residue in chain.get_residues():
                        if residue.get_id()[1] in resSeq:
                            resname_by_chain += seq1(residue.get_resname())
            resnames.append(resname_by_chain)

        self.config["binder"]["mutating_resnames"] = resnames
        return resnames

    def __generate_iteration_ID__(
        self, input_path: Path
    ) -> Tuple[str, List[str], List[List[int]], List[str]]:
        pdb_path = input_path / (self.config["main"]["name"] + ".pdb")
        chainIDs = self.config["binder"]["mutating_chainID"]
        resSeqs = self.config["binder"]["mutating_resSeq"]
        resnames = self.__get_mutating_resname__(pdb_path, chainIDs, resSeqs)

        iter_name = "".join(
            [f"-{chainID}_{resname}" for chainID, resname in zip(chainIDs, resnames)]
        )

        # Drop the leading -:
        return iter_name[1:], chainIDs, resSeqs, resnames

    def update_all_mdps(self, dir_path: Path, mdp_names: Dict) -> None:
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
