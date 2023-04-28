==========================================
Tutorial: running a simple optimization
==========================================

To recap, locuaz iteratively performs mutations on the selected residues of the binder, followed by minimization and
MD simulation. The trajectory is then scored with various scoring functions, to assess whether the mutation
results in higher or lower affinity. If the mutation results in better affinity, the mutation is accepted,
and the process is repeated. If the mutation does not significantly improve the affinity, then the mutants
are discarded and a new set of mutants are generated, based on the original complex(es).

In this tutorial, we are going to optimize the 3 CDRs from a nanobody against a p53 protein. This *p53* has a
*TWIST box* where the *TWIST* protein attaches and promotes its degradation. To prevent this, we're going to optimize
a nanobody that will compete for the binding of the *TWIST box*. For more info, check the `reference`_.

Assuming we have a conda environment named *locuaz*, let's start by activating it:

.. code-block:: console

    mamba activate locuaz

.. figure:: ./resources/simple_complex.png
        :alt: twistp53-nanobody complex

        Figure 1: snapshot of one optimized complex. **p53** is the yellow one on the left. Attached to its TWIST-box
        it's the nanobody, with a green backbone and with its CDRs 1, 2 and 3 colored magenta, orange and gray,
        respectively.


Starting files
----------------
This tutorial will provide a step-by-step guidance to the use of the locuaz protocol to optimize a nanobody (VHH)
as a p53 binder (figure shown above). The bulk of it will concentrate on the writing of the configuration yaml file.
A more detailed explanation of the available options, can be found in the :ref:`YAML configuration file`.

The materials for this tutorial are located in `examples/simple_tutorial`_ folder:

1.  ``d11.pdb``: the PDB file of the pre-equilibrated complex.
2.  ``config_simple.yaml``: the configuration file; all options go in here.
3.  ``mdp`` directory: minimization, NVT and NPT *GROMACS* input files.

.. _reference: http://dx.doi.org/10.1016/j.ccr.2012.08.003
.. _examples/simple_tutorial: https://istitutoitalianotecnologia-my.sharepoint.com/personal/walter_rocchia_iit_it/_layouts/15/onedrive.aspx?ga=1&id=%2Fpersonal%2Fwalter%5Frocchia%5Fiit%5Fit%2FDocuments%2FExamples&view=0

The configuration file
-----------------------
Options don't go all bunched up together, there are classified in sections. We will review each one of them.

paths
^^^^^^
If you are running this example, you will have to change every field on this section, since these are system dependent.
 * ``gmxrc``: is the location of *GROMACS* ``GMXRC`` script and binaries.
 * ``scoring_functions``: all scoring functions have to be in this folder, with a folder for each one and all its
   necessary files inside. So, for example, if you are using *bluues*, *evoef2* and *gmxmmpbsa*, inside the
   ``scoring_functions`` directory you will need the ``bluues`` folder with the **bluues** binary inside, the ``evoef2``
   with the **evoef2** binary inside and the ``gmxmmpbsa`` folder with the **gmxmmpbsa** input script inside.
 * ``mutator``: the necessary files for the mutator have to be here. We are using the *dlpr* mutator in this example,
   so the *DLPacker* files have to be here: **charges.rtp**, **DLPacker_weights.h5** and **library.npz**.
 * ``mdp``: *GROMACS* mdp input files have to be here.
 * ``input``: the starting PDB goes here.
 * ``work``: protocol generated files go here.

main
^^^^^
We are going to leave most values to their defaults and only set the name of our run.

 * ``name``: the name of the run. The input PDB located in ``[paths][input]`` has to be named after this name
   (``{name}.pdb``).
 * ``starting_epoch``: we are going to start on the default value of ``0``, but if your protocol run is a continuation
   of a previous one, you can set this value to other number in order to facilitate your posterior analysis.

protocol
^^^^^^^^
Global options of the run go here.

 * ``epochs``: the number of *epochs* we want to run. Remember that a failed *epoch*, that is, an *epoch* that fails
   to generate at least 1 *iteration* that improves the binding score is backed up (its folder is prefixed with ``bu_``)
   and is not included in the total number. So this will be the total number of successful epochs.
 * ``branches``: the number of parallel runs. If we look at the workflow from :ref:`Main idea`, it would be the 'width'
   of the workflow. See below for more info.
 * ``memory_size``: we want to prevent *locuaz* from mutating a position that was recently mutated, so we set this
   number to ``4``, this means that if position, say, ``128`` is mutated on epoch ``12``, then it won't be mutated again
   at least until epoch ``17``.
 * ``failed_memory_size``: similar to ``memory_size`` but it's only filled if the mutation of the position failed to
   improve affinity. Useful when we don't want the protocol to go back to a failed position for a long time, but at the
   same time we don't want to increase the ``memory_size`` too much, which would eliminate a lot of randomness from out
   run. We will set it to ``8``.

``constant_width`` is defaulted to ``True``, this means that each *epoch* will have 4 *iterations* so if, for example,
1 *complex* moves on to the next *epoch*, then 4 mutations will be performed on this complex, but if it were 3 complexes
then 2 of them would be mutated just once, and only one of them (chosen randomly), will be mutated twice; thus giving 4
new *iterations*.
If ``constant_width`` was ``False``, then ``branches`` is the number of mutations performed on each complex from the
previous step. Eg: if 2 complexes move on to the next epoch and ``branches=4``, then the next *epoch* will run 8
iterations, since 4 new complexes were obtained from each surviving complex.

generation
^^^^^^^^^^^
Now we begin to deal with a *locuaz* concept, :ref:`Units`. These are the moving parts of *locuaz*. The first one is the
mutation generator, the *unit* that is in charge of taking the sequence of the current complex and generating a new
sequence from it.

 * ``generator``: we are using the :ref:`SPM4gmxmmpbsa` generator, so later we will have to include *gmxmmpbsa* as a
   scoring function, so this generator can read the energy decomposition file from *gmxmmpbsa* and choose the position
   with the lowest contribution to the affinity as the position to mutate.
 * ``probe_radius``: this parameter is only used when the generator includes interface information, which is the case
   for SPM4gmxmmpbsa and others (eg: :ref:`SPM4i`). The *generator* uses *freesasa* to determine the CDR residues
   that form part of the interface and only considers those as potentials candidates for mutation. Since *freesasa* is
   a rolling-probe method, ``probe_radius`` allows the user to set the size of this probe. In this example we are using
   a radius of ``3``, a rather large probe, so more residues end up being classified as part of the interface.

Check :ref:`Mutation Generators` for a reference of the implementation, and :ref:`config['generation']` in the
:ref:`YAML configuration file` page for more details.

mutation
^^^^^^^^
This is another *unit*, the one that is in charge of performing the actual mutation.

 * ``mutator``: the external program to mutate the complex and find a suitable side-chain orientation. We are using
   ``dlpr`` since depends on the *DLPacker* program, which comes built-in with *locuaz* and also performs a nice
   reconstruction of the surrounding side-chains.
 * ``reconstruct_radius``: residues below this distance from the mutated position will also get their side-chains
   reoriented.

Check :ref:`Mutators` for a reference of the implementation, and :ref:`config['mutation']` in the
:ref:`YAML configuration file` page for more details.

pruning
^^^^^^^^
In this *unit*, you can set how the top *iterations* from an *epoch* will be selected to pass onto the next one.

 * ``pruner``: the external program to mutate the complex and find a suitable side-chain orientation. We are using
   ``dlpr``

md
^^^^
Yet anothe

target
^^^^^^^^
binder
^^^^^^^^
^^^^^^^^
scoring
^^^^^^^^

