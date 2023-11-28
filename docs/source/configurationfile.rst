========================
YAML configuration file
========================

locuaz behaviour is highly customizable, there are a lot of moving parts and alternative ways of carrying out the
same function. Hence this reference, in which we'll refer to the input configuration as ``config`` and explain each
of its "sections".

``config`` sections
--------------------
Given the high number of options, ``config`` has a hierarchical structure, where similar options are
grouped in sections.

``config['paths']``
^^^^^^^^^^^^^^^^^^^^^^^^
All paths go in this section. And all of them are required, except for:

1. ``config['paths']['input']``: only used when starting a new protocol run, since during restarts the protocol
   will read the necessary files from the current epoch.
   The input PDB file has to be here and be named ``{config['paths']['name']}.pdb``
2. ``config['paths']['tleap']``: only necessary when using tleap to build topologies. The **tleap** script and
   any other auxiliary file (``.frcmod`` or ``.lib`` files) should be here.

Some other things to highlight:

1. ``config['paths']['gmxrc']``: path to where the ``GMXRC`` is located, along with the GROMACS binary, usually
   called `gmx`
2. ``config['paths']['scorers']``: root directory where each scorer will have its own folder.
   Check :ref:`scorers:scorers` for more info.
3. ``config['paths']['mutator']``: mutator binary and/or parameters have to be here and be named appropriately.
   Check :ref:`mutators:Mutators` for more info.
4. ``config['paths']['work']``: name of the working dir. If it's an existing directory, the protocol will assume
   it's restarting from a previous run, if not, it will start a new one.


``config['main']``
^^^^^^^^^^^^^^^^^^^^^^^^
1. ``config['main']['name']``: system's name.
2. ``config['main']['mode']``: in case you want to just finish the MD run of the last epoch on your work dir, set
   this to ``run``, if you just want to score it set it to ``score``. In most scenarios, the default option ``evolve``
   is what you want. This is also the only ``config`` option that may be overrided by a CLI option
   (eg: ``--mode evolve``).
3. ``config['main']['prefix']``: the prefix files from the NPT run will get. Not very useful.
4. ``config['main']['starting_epoch']``: useful when a starting a new protocol from a PDB, or set of PDBs, that
   was already optimized by a protocol run, eg: when running an unrestrained optimization after a restrained one.

``config['protocol']``
^^^^^^^^^^^^^^^^^^^^^^^^

1. ``config['protocol']['epochs']``: .
2. ``config['protocol']['branches']``: .
3. ``config['protocol']['prevent_fewer_branches']``: .
4. ``config['protocol']['memory_size']``: .
5. ``config['protocol']['memory_positions']``: .
6. ``config['protocol']['failed_memory_size']``: .
7. ``config['protocol']['failed_memory_positions']``: .
8. ``config['protocol']['memory_aminoacids']``: not yet implemented.

``config['generation']`` (deprecated)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
1. ``config['generation']['generator']``: check the currently available generators below.
2. ``config['generation']['probe_radius']``: some generators exclude positions that are not currently on the
   interface, which is calculated by the rolling probe method
   (`freesasa <https://freesasa.github.io/doxygen/Geometry.html>`_); the bigger the probe, the more residues that
   will be classified as being part of the interface.

.. _config_creation:

``config['creation']``
^^^^^^^^^^^^^^^^^^^^^^
**Site selection options**

1. ``config['creation']['sites']``: number of sites to mutate. Each new branch
   will still get 1 mutation, so increasing this number will increase the number
   of new branches, if ``constant_width=false``.
2. ``config['creation']['sites_interfacing']``: only consider sites that are on
   the interface with the target.
3. ``config['creation']['sites_interfacing_probe_radius']``: a higher probe radius
   will increase the number of residues that are considered as being part of the
   interface.
4. ``config['creation']['sites_probability']``: set it to ``uniform`` if you want
   all positions to have the same chance of being chosen, or set it to ``mmpbsa``
   if you're already using the *mmpbsa* scorer and want to choose the sites that
   contribute the less to the interaction.
   Remember that to do this you need to add a section for residue decomposition
   to the ``gmxmmpbsa`` script file::

    /
    &decomp
    idecomp=2, dec_verbose=0,
    print_res="within 4"
    /

Check the :ref:`scorers:gmxmmpbsa` section for more info.

**Amino acid selection options**

These options affect the probability of each amino acid being chosen to be placed
at the already selected site.

1. ``config['creation']['aa_bins']``: list of strings where each element is a bin,
   represented as string of consecutive one-letter coded amino acids.
   If you don't want to group amino acids, just set it to 1 bin with all amino acids
   like this: ``[CDESTAGIMLVPFWYRNQHK]``
2. ``config['creation']['aa_bins_criteria']``: set it to ``without`` so amino acids
   will be chosen from all other bins, but the one that contains the current amino
   acid, before the mutation. The opposite effect is obtained when set to ``within``,
   only amino acids contained in the same bin as the current one can be chosen.
3. ``config['creation']['aa_probability']``: it can be set to either ``uniform``,
   ``ReisBarletta``, to use the probabilities extracted from the Reis & Barletta et. al.
   paper, and ``custom``, to set your own. In this last case, you'll also have to
   set the following option.
4. ``config['creation']['aa_probability_custom']``: a 20-element dictionary with
   the probability assigned to each amino acid.

``config['mutation']``
^^^^^^^^^^^^^^^^^^^^^^^^
1. ``config['mutation']['mutator']``: check the currently available mutators below.
2. ``config['mutation']['reconstruct_radius']``: when using the **dlpr** mutator, residues within this radius from
   the mutated residue will get their sidechains reoriented by `DLPacker <https://github.com/nekitmm/DLPacker>`_.
3. ``config['mutation']['allowed_nonstandard_residues']``: if there're non-protein residues
   and these are on the optimized interface, they need to be taken into account when
   optimizing the sidechain of the newly mutated residue. Add their resnames here.

``config['pruning']``
^^^^^^^^^^^^^^^^^^^^^^^^
1. ``config['pruning']['pruner']``: check the currently available pruners below.
2. ``config['pruning']['remaining_branches']``: you can set this value when the chosen pruner leaves a fixed number
   of branches after pruning.


``config['md']``
^^^^^^^^^^^^^^^^^^^^^^^^

``config['target']``
^^^^^^^^^^^^^^^^^^^^^^^^

``config['binder']``
^^^^^^^^^^^^^^^^^^^^^^^^

``config['scoring']``
^^^^^^^^^^^^^^^^^^^^^^^^

``config['statistics']``
^^^^^^^^^^^^^^^^^^^^^^^^

schema.yaml
------------

Input configuration files are validated against the following schema, which also works as a reference you can check
when in doubt, given its plain-english syntax. For example, you can check whether an option is mandatory or not
(``required``), if it requires a string, a number, etc. (``type``), if it has a default value (``default``), etc.

.. literalinclude :: ../../locuaz/schema.yaml
   :language: yaml
