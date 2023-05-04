===================================
Tutorial: using Tleap topologies
===================================

While the protocol uses *GROMACS* to perform MD simulation. It can also use **ambertools** to build an amber topology,
which can then be converted into *GROMACS* topology. This allows the use of force fields that are not available in *GROMACS*,
the inclusion of ligands, non-standard amino acids, etc.

In this example, we are going to optimize a nanobody towards a protein that contains Zinc and coordinates it with
amino acids that can't be represented on a regular force-field. Hence, we're going to need ***Tleap*** to build the
topology and **parmed** to turn it into something *GROMACS* can work with. Both these tools come with **ambertools**,
which comes with the protocol. For mor info, check the :ref:`installation:Installation` section.

.. figure:: ./resources/tleap_complex.png
        :alt: p53-nanobody complex

        Figure 1: snapshot of one optimized complex. **p53** is the yellow one on the left, with its loops colored red and
        blue, these loops have to be stabilized so it doesn't loose its function; the zinc atom and its coordinating
        residues are on the bottom-left corner. The nanobody is the green one on the right, with its CDRs 1, 2 and 3
        colored magenta, orange and gray, respectively.

As always, the name of our environment is *locuaz*, so we start by activating it.

.. code-block:: console

    mamba activate locuaz

Setting up the system
----------------------
The starting complex (**p53-VHH**) in this tutorial has been obtained using `HADDOCK`_, followed by an
equilibration. In this example, the 4-coordinated zinc metal center in **p53** is described with Zinc
AMBER force field **ZAFF**. Therefore, the topology is built with a ***Tleap*** script, which requires
additional parameter files ``ZAFF.frcmod`` and ``ZAFF.prep``.

The topology of the complex provided in this example has been prepared with the following *Tleap* script:

.. code-block:: console

    source oldff/leaprc.ff99SBildn
    source leaprc.water.tip3p
    addAtomTypes { { "ZN" "Zn" "sp3" } { "S2" "S" "sp3" } { "N1" "N" "sp3" } }
    loadamberparams frcmod.ions1lm_126_tip3p
    loadamberprep ZAFF.prep
    loadamberparams ZAFF.frcmod
    mol = loadpdb nb.pdb

    bond mol.195.ZN mol.81.SG
    bond mol.195.ZN mol.84.ND1
    bond mol.195.ZN mol.143.SG
    bond mol.195.ZN mol.147.SG
    bond mol.217.SG mol.292.SG     ##adding the S-S bond

    solvatebox mol TIP3PBOX 10.0
    addions mol CL 0
    addions mol NA 0

    saveamberparm mol amb_nb.prmtop amb_nb.inpcrd
    quit

The script can be run as:

.. code-block:: console

    tleap -f tleap

The topology was then converted into the *GROMACS* topology format using `ParmEd`_, `acpype`_ is
an alternative but we recommend staying with *ParmEd* since it's the same *locuaz* uses internally.

Lastly, the system was minimized and equilibrated using the same *mdp* *GROMACS* input files we will
be using during the protocol, ``min.mdp`` and ``nvt.mdp``. For the *NPT* run, we ran for 10ns and saw
no changes in the interface. Ideally you do not want to start with a complex that changes its interface
too much, since this will change the affinity and the scoring of the mutations loose meaning.

Necessary files
----------------
As in :ref:`tutorialsimple:Tutorial: running a simple optimization`, we're going to focus on the writing
of the YAML config file. A more detailed explanation of the available options, can be found in the
:ref:`configurationfile:YAML configuration file`. The materials for this tutorial are located in
the ``examples/tleap_tutorial`` folder:

1. ``nb.pdb``: the PDB file of the pre-equilibrated complex.
2. ``tleap``: *Tleap* script to build the topology of the system each time a mutation is performed. This
   script will be identical to the one above, with the exception of the ``solvatebox`` line, since the
   solvent is already there. Another thing to notice is the usage of ``addions``. We keep this commands
   since *Tleap* will be responsible of keeping neutrality of the system and avoid using ``addions2`` since
   we need it to replace water molecules each time it ads ions, to keep the *N* of the system constant.
3. ``ZAFF.frcmod`` and ``ZAFF.prep`` (auxiliary Zn parameters)
4. ``config_tleap.yaml``: the input file to run the protocol.
5. ``mdp`` directory: minimization, NVT and NPT *GROMACS* input files.


The configuration file
-----------------------
We will focus on the new options that didn't show up on :ref:`tutorialsimple:Tutorial: running a simple optimization`.

paths
^^^^^^
.. code-block:: console

    paths:
        gmxrc: /apps/*GROMACS*/2021.4/gcc7-ompi4.1.1-cuda11.1-plm2.8.0/bin
        scoring_functions: /work/rtandiana/mdp/SF
        mutator: /work/rtandiana/mdp/SF/dlpacker
        *Tleap*: /work/rtandiana/Optimization/New-ZAFF/NB112/C9/input
        mdp: /work/rtandiana/mdp
        input: [ /work/rtandiana/Optimization/New-ZAFF/NB112/C9 ]
        work:  /work/rtandiana/Optimization/New-ZAFF/NB112/C9/work_dir

 * *Tleap*: the path to the ***Tleap*** scripts. It is mandatory if *Tleap* is used.

main
^^^^^

In the main sections, the name of the PDB files are defined, and it has to match the pdb file provided in the input directory. The running mode of the protocol is set to evolve.

.. code-block:: console

    main:
        name: nb
        mode: evolve

protocol
^^^^^^^^
In the protocol section, several important options concerning the protocol have to be specified.

 * epochs: The number of epochs desired
 * branches: The number of iterations at each epochs, which usually correlate to the number of GPUs available. Each iteration corresponds to different target mutation
 * prunner: The method adopted to pick the best iteration(s) in each epoch
 * generator: The algorithm to generate the mutation
 * mutator: The algorithm to generate the mutated structure
 * memory_size: The number of selected position of mutation of previous epochs that the protocol will retain
 * failed_memory_size: The number of selected position of mutation of previous failed epochs that the protocol will retain

The memory_size and failed_memory_size options assist to prevent the protocol to perform mutation at previously mutated residues.

.. code-block:: console

    protocol:
        epochs: 20
        branches: 4
        prunner: threshold
        generator: SPM4i
        mutator: dlpr
        memory_size: 4
        failed_memory_size: 4

generation
^^^^^^^^^^^

mutation
^^^^^^^^


pruning
^^^^^^^^


md
^^^^
.. code-block:: console

    md:
        gmx_mdrun: gmx mdrun
        mdp_names:
            min_mdp: min.mdp
            nvt_mdp: nvt.mdp
            npt_mdp: npt.mdp
        ngpus: 4
        mpi_procs: 1
        omp_procs: 8
        pinoffsets: [ 0, 32, 64, 96 ]
        use_tleap: true

 * gmx_mdrun: The *GROMACS* command to perform mdrun
 * mdp_names: The name of the mdp files present in the mdp folders specified above
 * ngpus: The number of GPUs available
 * mpi_procs: typically 1
 * omp_procs: The number of threads used for each MD runs
 * pinoffsets: Pinning the threads to specific positions to maximize the performance. This values depend on the GPU architecture
 * use_tleap: True, this option is specified only if *Tleap* is used to build the topology.

target
^^^^^^^^
.. code-block:: console

    target:
        chainID: [A]

binder
^^^^^^^^
.. code-block:: console

    binder:
        chainID: [B]
        mutating_chainID: [B,B,B]
        mutating_resSeq: [[220,221,222,223,224,225,226,227],[248,249,250,251,252,253,254],[294, 295, 296, 297, 298, 299, 300]]
        mutating_resname: [[S,G,F,D,F,S,D,A],[R,S,G,L,A,T,S],[K,S,R,R,G,Q,G]]

In the binder section, where the single point mutation will be performed, the following options have to be specified:
    *	chainID: The chain IDs of the binder
    *	mutating_chainID: The chain IDs where the mutation is desired. Since the mutating sequence is listed separately for each CDRs, the chain ID has to be a list also.
    *	mutating_resSeq: The residue sequence of the desired mutation sites. In this example, the sequences are listed separately for each CDRs.
    *	mutating_resname: The amino acid residues in one letter format, that correspond to the mutating_resSeq


scoring
^^^^^^^^
.. code-block:: console

    scoring:
        functions: [ bluuesbmf, piepisa, evoef2, gmx_mmpbsa ]
        consensus_threshold: 3
        nthreads: 80
        mpiprocs: 2
        start: 50
        end: -1

2 new options show up with respect to the :ref:`tutorialsimple:Tutorial: running a simple optimization`

 * ``start``: Useful if you want to skip a few frames before starting to score. 0-indexed.
 * ``end``: Also 0-indexed. Defaults to ``-1``, which means all remaining frames.


Running the protocol
---------------------
There's nothing new here with respect to the simple tutorial, we just run the protocol with our config file

.. code-block:: console

    mamba activate locuaz
    python /home/user/locuaz/locuaz/protocol.py config_tleap.yaml


And as always, the protocol will create the working directory folder and inside of it, a folder for
each *iteration*:

.. figure:: ./resources/tleap_workdir.png
        :alt: directory structure of an iteration folder

        Figure 2: the look of any *iteration* folder after it has been finished. *Tleap* related files
        are highlighted.


.. _HADDOCK: https://wenmr.science.uu.nl/haddock2.4/
.. _ParmEd: https://github.com/ParmEd/ParmEd
.. _acpype: https://github.com/alanwilter/acpype