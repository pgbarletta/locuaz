from pathlib import Path
from projectutils import WorkProject, Epoch, Iteration
from fileutils import FileHandle, DirHandle
from biobb_model.model.mutate import mutate
from molecules import GROComplex
from mutator import generate_new_binders
from gromacsutils import remove_overlapping_waters

# This is a specific protocol, others will be added
def initialize_new_epoch(work_pjct: WorkProject):
    old_epoch = work_pjct.epochs[-1]
    current_epoch = Epoch(old_epoch.id + 1, {})

    for old_iter_name, old_iter in old_epoch.items():

        (
            mut_resSeq,
            mut_strings,
            new_iterations_names,
            new_iterations_resnames,
        ) = generate_new_binders(old_iter, width=work_pjct.config["main"]["width"])

        for mut_string, iter_name, iter_resnames in zip(
            mut_strings, new_iterations_names, new_iterations_resnames
        ):
            iter_path = Path(
                str(work_pjct.dir_handle), Path(str(old_epoch.id + 1) + "-" + iter_name)
            )
            this_iter = Iteration(
                DirHandle(iter_path, make=True),
                iter_name=iter_name,
                chainIDs=old_iter.chainIDs,
                resnames=iter_resnames,
                resSeqs=old_iter.resSeqs,
            )
            mut_name = "mut_" + work_pjct.config["main"]["name"]
            mut_pdb = iter_path / (mut_name + ".pdb")
            props = {"mutation_list": mut_string}
            mutate(
                input_pdb_path=str(old_iter.complex.pdb.file.path),
                output_pdb_path=str(mut_pdb),
                properties=props,
            )

            # Create complex with coordinates (from npt run) and topology (should be zip)
            # This mutated PDB lacks hydrogens in the mutated residue, so it won't
            # have the same topology as the gro and top files.
            overlapped_cpx = GROComplex.from_pdb(
                name=mut_name,
                input_dir=iter_path,
                target_chains=work_pjct.config["target"]["chainID"],
                binder_chains=work_pjct.config["binder"]["chainID"],
                gmx_bin=work_pjct.config["md"]["gmx_bin"],
            )

            this_iter.complex = remove_overlapping_waters(
                work_pjct.config, overlapped_cpx, mut_resSeq
            )

            current_epoch[iter_name] = this_iter
    work_pjct.new_epoch(current_epoch)
