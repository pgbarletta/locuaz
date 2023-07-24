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
2. ``config['paths']['scoring_functions']``: root directory where each scorer will have its own folder.
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

``config['generation']``
^^^^^^^^^^^^^^^^^^^^^^^^
1. ``config['generation']['generator']``: check the currently available generators below.
2. ``config['generation']['probe_radius']``: some generators exclude positions that are not currently on the
   interface, which is calculated by the rolling probe method
   (`freesasa <https://freesasa.github.io/doxygen/Geometry.html>`_); the bigger the probe, the more residues that
   will be classified as being part of the interface.


``config['mutation']``
^^^^^^^^^^^^^^^^^^^^^^^^
1. ``config['mutation']['mutator']``: check the currently available mutators below.
2. ``config['mutation']['reconstruct_radius']``: when using the **dlpr** mutator, residues within this radius from
   the mutated residue will get their sidechains reoriented by `DLPacker <https://github.com/nekitmm/DLPacker>`_.

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
