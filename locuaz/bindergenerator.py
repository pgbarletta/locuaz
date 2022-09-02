from typing import List, Tuple
from projectutils import Iteration
from random import choice, sample
from mutator import Mutation
import logging


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


def generate_mutations(iteration: Iteration, *, branches: int) -> List[Mutation]:
    if branches > 19:
        logging.warning(
            f"{branches} is over 19 but this mutator generates "
            "mutations for 1 position, hence, 19 mutations is the maximum number "
            "of possible mutations. "
        )
        branches = 19
    # First, choose the position to mutate:
    n_chains = len(iteration.chainIDs)
    idx_chain = choice(range(0, n_chains))
    n_residues = len(iteration.resSeqs[idx_chain])
    idx_residue = choice(range(0, n_residues))
    # Then, save the current AA
    old_aa = iteration.resnames[idx_chain][idx_residue]

    n_cat = len(CAT_AA1_LIST)
    categories = set(range(n_cat))
    new_aas = {old_aa}
    while branches != 0:
        cat_idx = choice(tuple(categories))
        categories.difference_update({cat_idx})

        shuffled_aas = sample(CAT_AA1_LIST[cat_idx], len(CAT_AA1_LIST[cat_idx]))
        for new_aa in shuffled_aas:
            if new_aa not in new_aas:
                new_aas.add(new_aa)
                branches -= 1
                break
        else:
            raise RuntimeError(
                "Can't generate new binder. This is a logic error. "
                "This shouldn't happen."
            )

        if len(categories) == 0:
            # All categories have already been chosen from `N` times.
            # Allow all of them again for the `N+1` iteration, except those
            # that are exhausted already
            categories = set(range(n_cat))
            for i in range(n_cat):
                if set(CAT_AA1_LIST[i]).issubset(new_aas):
                    # All AAs from thise category have already been chosen
                    categories.difference_update({i})

    # Finally, build the list of mutation objects
    new_aas.difference_update({old_aa})
    mut_chainID = iteration.chainIDs[idx_chain]
    mut_resSeq = iteration.resSeqs[idx_chain][idx_residue]
    mutations = [
        Mutation(
            chainID=mut_chainID,
            resSeq=mut_resSeq,
            old_aa=old_aa,
            new_aa=aa,
            chainID_idx=idx_chain,
            resSeq_idx=idx_residue,
        )
        for aa in new_aas
    ]

    return mutations
