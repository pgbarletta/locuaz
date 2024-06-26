=======
History
=======

0.7.5.3 (2024-1-)
------------------
 * Fix install. Add networkx missing dependency.

0.7.5.2 (2024-1-)
-------------------
 * Fix wrong epoch id on branches when using roundrobin pruner.
 * Fix ``get_interface_surface()`` when 'CYX' residues are present. When
   collecting resnames from freesasa, use AA_MAP to turn 'CYX' into 'CYS'.
 * Fix ``get_pdb_tpr()`` when an ion is assigned its own chainID.
   'gmx trjconv', instead of reading the chainIDs from the .tpr file, generates them
   by itself, leaving the ion without its chainID and misnaming the chainIDs of
   the following molecules. Added a new method ``get_pdb()`` that can use biobb
   or MDAnalysis to build the PDB. This will be a temporary fix until GROMACS
   changes its behaviour, or a permanent one (much more likely).
 * Fix ``schema.yaml`` to allow for resnames shorter than 3 letters when setting
   ``allowed_nonstandard_residues``.
 * Fix ``fixup_top()`` better error reporting.

0.7.5.1 (2023-12-)
------------------
 * epoch initialization: fix broken initialization when tleap wasn't being used
   by adding locks around pdb2gmx call, which apparently can't be parallelized.

0.7.5 (2023-12-)
------------------
 * MD runs: add tolerance to failed MDs and better login when a run fails.
 * MD runs: only run NVT after all MIN are done. MD issues were coming up due to
   the asynchronicity of these 2 steps.
 * epoch initialization: make branch creation parallel.

0.7.4.8 (2023-12-)
------------------
 * scoring: silence warning from writing fixed PDB in image_traj()
 * MD runs: fix run_min_nvt_epoch(). It was serializing NVT runs instead of parallelizing them.
 * Mutator: fix evoef2 initialization
 * Mutation: warn if setting ``allowed_nonstandard_residues`` and using any other
   Mutator than dlpr, since ``dlp`` may create overlaps and ``evoef2`` doesn't see ligands.

0.7.4.7 (2023-12-)
------------------
 * Mutation Creation: fix ``uniform.csv`` for uniform amino acid selector
   probability.
 * Mutation Creation: don't accept amino acids with zero probability.
 * branch creation: fix remove_overlapping_solvent's behaviour when no water is added.

0.7.4.6 (2023-12-)
------------------
 * Fix ``uniform.csv`` for uniform amino acid selector probability.
 * Fix epoch initialization. If the branch creation loop fails twice, it will
   abort branch creation and just continue with the protocol.

0.7.4.5 (2023-12-)
------------------
 * Fix auto_md_params() again.
 * Fix bug that would come up during an epoch initialization when top branches
   only differ by 1 amino acid, that is in the same position that is currently
   being mutated. If the mutation creator assigns the same amino acid to both of
   them, identical branches would be generated.
   This scenario is rare, but may come up during the first steps of an optimization
   that is using little to no memory of past mutation positions and hence may
   mutate the same position twice in a row or with few epochs in between.
   *locuaz* allows the creator repeat amino acids for different branches since it
   assumes that branches will differ in other positions.
   To remediate this bug, *locuaz* will generate new mutations from other random
   top branches to replace those branches that were repetitions.

0.7.4.4 (2023-11-)
------------------
 * Remove ``OMP_NUM_THREADS`` environment variable.
 * Fix auto_md_params().

0.7.4.3 (2023-11-)
------------------
 * Force set omp threads environment variable at each MD run.

0.7.4.2 (2023-11-)
------------------
 * Fix NUMA regions.

0.7.4.1 (2023-11-)
------------------
 * Allow more NUMA regions than GPUs available.

0.7.4 (2023-11-)
------------------
 * Fix bug in creation options schema.
 * Fix bug ``BaseMutator.port_mutation()`` when a non standard residue is surrounding
   the mutated residue.
 * Fix bug when host sets ``OMP_NUM_THREADS`` environment variable and ``mps=true``.

0.7.3 (2023-11-)
------------------
 * Add ``allowed_nonstandard_residues`` to ``mutation`` options.

0.7.2 (2023-11-)
------------------
 * Fix broken install due to pygraphviz being broken on pypi.

0.7.1 (2023-11-)
------------------
 * Allow mutating more than 1 site.
 * Add pygraphviz dependency to graph the DAG.
 * Fix input config and PDB checking when converting 3-letter code AAs into
   1-letter code.
 * Fix error message when converting 3-letter code AAs into 1-letter code.
 * Started adding reformat with ruff.

0.7.0 (2023-10-)
------------------
 * Added warning when 'autodockvina' scorer is used what no resname was set in
   'allowed_nonstandard_residues'. The former is usually used to score
   interactions with small molecules that will be discarded from the PDBs used
   for scoring, unless their resnames show up on the
   'allowed_nonstandard_residues' list.
 * Added RoundRobin pruner. It'll take the current branches and the top branches
   from the previous epoch and select ``N`` branches as the new top branches.
   As a consequence, failed epochs won't be branded as such and branches from an
   epoch ``i`` may come from a mix of branches from the epochs ``i-1`` and
   ``i-2``.
 * Added MutationCreator as a future replacement of MutationGenerator. Favouring
   composition over inheritance, MutationCreator is fully user-customizable
   instead of offering a set of fixed options as MutationGenerator.
   MutationCreator offers all the possibilities from MutationGenerator and more.

0.6.3 (2023-10-)
------------------
 * Rename positional restraints from "target", "binder" and "rest" to "posres"
   and "posres_water".

0.6.2 (2023-09-)
------------------
 * Pin gmx-mmpbsa to 1.6.1 since 1.6.2 pins pandas to 1.2.2 which is broken.

0.6.1 (2023-09-)
------------------
 * Support Python version 3.10 and onwards.

0.6.0 (2023-09-)
------------------
 * Fixed bug when NPT positional restraints weren't used.
 * Support Python version 3.11 and onwards.

0.5.3 (2023-07-)
------------------
 * Renamed ``scoring functions`` to ``scorers``.
 * Added support for positional restraints.
 * Pinned Python version to 3.10.X.

0.5.2 (2023-07-)
------------------
 * Fix the PDB left as  reference inside the ``scoring`` dir. ``fix_npt_{name}.pdb`` is left as a topology
   for the cleaned trajectory file ``fix_npt_{name}.xtc``. Now it contains PDB contains **chainID** info.
 * Pinned Python version to 3.10 and newer.

0.5.1 (2023-07-)
------------------
 * Fix *DLPacker* data download through pip.

0.5.0 (2023-07-)
------------------
 * Added MPS usage. Now multiple runs can be queued up onto the same GPU and *locuaz* will decide the parameters for
   each process (which GPU to use, how many threads for OMP and for MPI and the pinoffset for the run).
   Expected improved throughput: ``1.3-2.0``.
 * Added support for positional restraints when building topology with *tleap* by defining ``-DPOSRES_TARGET``
   to restrain the target, ``-DPOSRES_BINDER`` to restrain the binder and ``-DPOSRES`` for everything else.
 * Removed ``prefix`` option to set a custom prefix to the files generated by the NPT run.
   Now the prefix is always ``"npt_"``
 * Added resiliency against uninitialized current epoch. If one of the current branches doesn't have the initial PDB,
    GRO, ZIP and TPR files, then the whole epoch is backed up on cli.py and the protocol will later initialize a
    whole new epoch.
 * Fixed ``gmxmmpbsa`` bug when MPI was not used.
 * Better plot for the DAGS at ``graphs.png``
 * Better login.
 * Added **Developing** section to the reference docs.

0.4.1 (2023-06-)
------------------
 * Renamed ``Iteration`` abstraction to ``Branch``
 * Made ``previous_branches``, ``current_branches``, ``top_branches`` variables in the tracking file ``tracking.pkl``
   relative paths to the work dir. This allows the work dir to be moved around without errors.

0.3.9 (2023-06-)
------------------
 * Added ``locuaz`` as executable.

0.3.8 (2023-05-)
------------------
 * *DLPacker* data files ``library.npz`` and ``charges.rtp`` are now inscluded with the install. Only the weights have
   to be downloaded and extracted into a dir whose path must be specified in the ``config['paths']['mutator']`` option.

0.3.7 (2023-05-12)
------------------
 * Added Directed Acyclic Graph tracking of the protocol, so a plot of the progression of the protocol can be done,
   both of the branch names and the mutations performed on each mutation.
 * Added docs on https://locuaz.readthedocs.io/
 * Made DLPacker part of the repo. Used for performing mutations.
 * Added metropolis Pruner.

0.2.1 (2023-04-20)
------------------
* The protocol is now fully installable by pip, provided that ambertools and tensorflow are present in the conda environment (no available pip install for them)

0.2.0 (2022-05-13)
------------------
* First fully functional release.

0.1.0 (2022-05-25)
------------------
* First release on PyPI.
