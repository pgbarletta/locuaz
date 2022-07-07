#!/usr/bin/env python

"""Main module."""
import cli
import sys
import os
from io import StringIO
from pathlib import Path
from fileutils import FileHandle
import projectutils as pu
from molecules import GROTopology, GROStructure, GROComplex, PDBStructure
from biobb_model.model.mutate import mutate
from gromacsutils import write_chain_selection_ndx, write_non_overlapping_ndx
from run import run_min_nvt_npt, prepare_for_mdrun
from scoring import score, initialize_scoring_folder


def main():

    config = cli.main()
    work_pjct = pu.WorkProject(config)

    # Get the baseline scores
    # last_epoch = work_pjct.epochs[-1]
    # this_iter_name = next(iter(last_epoch))
    # prepare_for_mdrun(work_pjct, this_iter_name)
    # run_min_nvt_npt(work_pjct, this_iter_name, config)
    initialize_scoring_folder(work_pjct, next(iter(work_pjct.epochs[0].keys())), config)
    # score(work_pjct, this_iter_name, config)

    print("ASBHDKKJABHDSCAJSBHDC")

    # mut_pdb = this_iter.complex.dir_handle.dir_path / (
    #     "mut_" + this_iter.complex.pdb.file.path.name
    # )

    # prope = {"mutation_list": "B:Pro3Glu", "use_modeller": False}
    # mutate(
    #     input_pdb_path=str(this_iter.complex.pdb.file.path),
    #     output_pdb_path=str(mut_pdb),
    #     properties=prope,
    # )

    # resid = 3
    # binder_ndx = write_chain_selection_ndx(
    #     "binder",
    #     ["B"],
    #     PDBStructure.from_path(mut_pdb),
    #     this_iter.dir_handle.dir_path,
    # )

    # ndx_handle, wat_cnt = write_non_overlapping_ndx(mut_pdb, binder_ndx, resid)

    # # gmx trjconv -f mut_1ppg.pdb -s mut_1ppg.pdb -o nonwat_mut_1ppg.pdb -n non_overlapping.ndx

    # print(wat_cnt)


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
