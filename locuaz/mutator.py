from pathlib import Path
from molecules import PDBStructure, GROStructure, ZipTopology, GROComplex
from typing import List, Tuple
from projectutils import Iteration
from biobb_model.model.mutate import mutate
from random import choice
from Bio.SeqUtils import seq1, seq3

# fmt: off
AA3_LIST = ("Ala", "Art", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "His", "Ile",
    "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp",  "Tyr", "Val")
AA1_LIST = ('A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
    'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')

NEG_AA1_LIST = ('D', 'E', 'S', 'T')
POS_AA1_LIST = ('R', 'N',  'Q',  'H', 'K')
PHO_AA1_LIST = ('A', 'I', 'M', 'F', 'L', 'W', 'Y', 'V')
MIS_AA1_LIST = ('C', 'G', 'P')
CAT_AA1_LIST = (NEG_AA1_LIST, POS_AA1_LIST, PHO_AA1_LIST, MIS_AA1_LIST)
# fmt: on


def generate_new_binders(
    iteration: Iteration, *, width: int
) -> Tuple[int, List[str], List[str], List[List[str]]]:
    # First, choose the position to mutate:
    n_chains = len(iteration.chainIDs)
    idx_chain = choice(range(0, n_chains))
    n_residues = len(iteration.resSeqs[idx_chain])
    # idx_residue = choice(range(0, n_residues))
    # TODO: remove this when trying out an appropiate system:
    idx_residue = choice(range(2, n_residues - 2))
    old_aa1 = iteration.resnames[idx_chain][idx_residue]
    # mutated resSeq:
    mut_resSeq = iteration.resSeqs[idx_chain][idx_residue]

    n_cat = len(CAT_AA1_LIST)
    new_aas1 = set()
    categories = list(range(n_cat))
    while width != 0:
        # Choose a category
        i = choice(categories)
        # Choose an AA from said category
        new_aa1 = choice(CAT_AA1_LIST[i])
        # Make sure it's actually new.
        if (new_aa1 != old_aa1) and (new_aa1 not in new_aas1):
            new_aas1.add(new_aa1)
            # Make sure this category of AA is not chosen again.
            categories = [cat for cat in categories if cat != i]
            if len(categories) == 0:
                # Renew the available categories to choose from when they
                # get exhausted
                categories = list(range(n_cat))
            width -= 1

    # Now, each mutation will result in a new iteration. For each one of them:
    ##### 1_ construct the string to perform the mutation,
    ##### 2_ generate the new iteration name,
    ##### 3_ generate the new resname, featuring the new amino acid.
    mut_texts = []
    new_iterations_names = []
    new_iterations_resnames = []
    for aa in new_aas1:
        mut_chainID = iteration.chainIDs[idx_chain]
        # Build the string that will go to biobb_model.model.mutate
        mut_texts.append(f"{mut_chainID}:{seq3(old_aa1)}{mut_resSeq}{seq3(aa)}")

        # Now, generate the new names for each iteration and resnames:
        iter_name = ""
        new_iteration_resnames = []
        for chainID, resname in zip(iteration.chainIDs, iteration.resnames):
            if chainID == mut_chainID:
                # This is the mutated chainID
                new_resname = resname[:idx_residue] + aa + resname[idx_residue + 1 :]

            else:
                # This one remains the same
                new_resname = "".join([residue for residue in resname])
            iter_name += f"-{chainID}_{new_resname}"
            new_iteration_resnames.append(new_resname)
        new_iterations_resnames.append(new_iteration_resnames)
        # Drop the leading '-'
        new_iterations_names.append(iter_name[1:])

    return mut_resSeq, mut_texts, new_iterations_names, new_iterations_resnames


# def build_grocomplex(
#     name: str, iter_path: Path, target_chains: List, binder_chains: List
# ) -> "GROComplex":
#     try:
#         str_pdb = PDBStructure.from_path(iter_path / ("npt_" + name + ".pdb"))
#         str_gro = GROStructure.from_path(iter_path / ("npt_" + name + ".gro"))
#         top = ZipTopology.from_path(iter_path / ("wet_topol.zip"))
#         top.update_chains(target_chains=target_chains, binder_chains=binder_chains)
#         # traj = XtcTrajectory.from_path(iter_path / ("npt_" + name + ".xtc"))
#         # tpr = TPRFile.from_path(iter_path / ("npt_" + name + ".tpr"))
#     except Exception as e:
#         print(
#             f"Could not get input files from: {iter_path}",
#             flush=True,
#         )
#         raise e
#     else:
#         cpx = GROComplex(name, iter_path, str_pdb, top, str_gro)
#         # cpx.tra = traj
#         # cpx.tpr = tpr
#         return cpx
