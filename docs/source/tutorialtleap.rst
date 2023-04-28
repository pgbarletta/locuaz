=================================
Tutorial: using tleap topologies
=================================

Introduction
------------

While the protocol uses *GROMACS* to perform MD simulation. It can also use **ambertools** to build an amber topology,
which can then be converted into *GROMACS* topology. This allows the use of force fields that are not available in *GROMACS*,
the inclusion of ligands, non-standard amino acids, etc.

In this example, we are going to optimize a nanobody towards a protein that contains Zinc and coordinates it with
amino acids that can't be represented on a regular force-field. Hence, we're going to need **tleap** to build the
topology and **parmed** to turn it into something *GROMACS* can work with. Both these tools come with **ambertools**,
which comes with the protocol. For mor info, check the :ref:`Installation` section.

.. figure:: ./resources/tleap_complex.png
        :alt: p53-nanobody complex

        Figure 1: snapshot of one optimized complex. **p53** is the yellow one on the left, with its loops colored red and
        blue, these loops have to be stabilized so it doesn't loose its function; the zinc atom and its coordinating
        residues are on the bottom-left corner. The nanobody is the green one on the right, with its CDRs 1, 2 and 3
        colored magenta, orange and gray, respectively.

As always, the name of our environment is *locuaz*, so we start by activating it.

.. code-block:: console

    mamba activate locuaz


Necessary files
----------------

As in :ref:`Tutorial: running a simple optimization`, we're going to focus on the writing of the YAML config file.
A more detailed explanation of the available options, can be found in the :ref:`YAML configuration file`.
The materials for this tutorial are downloaded along with the codes, within the ‘example’ folder, which are:

1.  nb.pdb (The PDB file of the pre-equilibrated complex)
2.  tleap (Tleap script to build the topology of the system, more details in the following section)
3.  ZAFF.frcmod and ZAFF.prep (auxiliary Zn parameters)
4.  :download:`config_tleap.yaml<../resources/config_tleap.yaml>` (The input file to run the protocol)
5.  mdp files (*GROMACS* parameter files)

Setting up the system
----------------------

In this example, we are using DLPacker as the mutator,

It is advisable to use a pre-equilibrated system as the starting complex for the optimization protocol.
In general, the following guidelines can be followed for the systems preparation:
1. The binder and target are docked with any docking program of choice, to generate multiple possible complex structures.
2. The selected or all complex structures are then equilibrated in a box of water with their MD program of choice.
3. The pdb file of the equilibrated structure can directly be used for the protocol.

The starting complex (p53/VHH) in this tutorial has been obtained by `
HADDOCK <https://wenmr.science.uu.nl/haddock2.4/>`_, followed by an equilibration.
In this example, the 4-coordinated zinc metal center in p53 is described with Zinc AMBER force field (ZAFF).
Therefore, the topology is built with a tleap script, which requires additional parameter files:
ZAFF.frcmod and ZAFF.prep.
The topology of the complex provided in this example has been prepared with the following tleap script:

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

Now, the topology has to be converted into the *GROMACS* topology format. Internally, locuaz uses
`ParmEd <https://github.com/ParmEd/ParmEd>`_ to do this, and we recommend to do the same.
Others may prefer to use `acpype <https://github.com/alanwilter/acpype>`_:

.. code-block:: console

    acpype -p amb_nb.prmtop -x amb_nb.inpcrd


Now, the minimization and 5ns of MD simulations can be performed with *GROMACS* to equilibrate the system, before continuing the optimization protocol.

Note that the protocol will maintain the size of the box given at the start. Therefore, in the tleap script provided to the protocol, the line "solvatebox mol TIP3PBOX 10.0" has to be removed. In addition, the addition of ions (either Na or Cl) at each iteration is taken care of by the tleap scripts.

Preparing the files
------------------------

The following files are needed to run the protocol, and their location should be specified in the input file (explained further later):
1.	tleap scripts and the additional parameter files (compulsory if tleap is used)
2.	The PDB file of the pre-equilibrated complex
3.	The input file for the protocol, with yaml extension. In this example, it is called config_tleap.yaml

In the input file, config_tleap.yaml, different options have to be specified:
1.	In the path sections, the paths to different folders have to be specified:
    *	gmxrc: the path to the *GROMACS* executable
    *	scoring_functions: folders containing the executable of different scoring functions (more details refer to the github page)
    *	mutator: folders containing the executable to generate mutated structures. In this example, DLPacker is used.
    *	tleap: the path to the Tleap scripts. It is mandatory if tleap is used.
    *	mdp: folder containing the *GROMACS* parameters
    *	input: folder containing the pdb files. Note that multiple files can be introduced as the starting structures, but in this example, we are using only 1 starting structure.
    *	work: The path where the working directory folder will be created, and where the results will be located. If it’s a new run, this directory should not exist.


.. code-block:: console

    paths:
        gmxrc: /apps/*GROMACS*/2021.4/gcc7-ompi4.1.1-cuda11.1-plm2.8.0/bin
        scoring_functions: /work/rtandiana/mdp/SF
        mutator: /work/rtandiana/mdp/SF/dlpacker
        tleap: /work/rtandiana/Optimization/New-ZAFF/NB112/C9/input
        mdp: /work/rtandiana/mdp
        input: [ /work/rtandiana/Optimization/New-ZAFF/NB112/C9 ]
        work:  /work/rtandiana/Optimization/New-ZAFF/NB112/C9/work_dir

2.	In the main sections, the name of the PDB files are defined, and it has to match the pdb file provided in the input directory. The running mode of the protocol is set to evolve.

.. code-block:: console

    main:
        name: nb
        mode: evolve

4.	In the protocol section, several important options concerning the protocol have to be specified.
    *	epochs: The number of epochs desired
    *	branches: The number of iterations at each epochs, which usually correlate to the number of GPUs available. Each iteration corresponds to different target mutation
    *	prunner: The method adopted to pick the best iteration(s) in each epoch
    *	generator: The algorithm to generate the mutation
    *	mutator: The algorithm to generate the mutated structure
    *	memory_size: The number of selected position of mutation of previous epochs that the protocol will retain
    *	failed_memory_size: The number of selected position of mutation of previous failed epochs that the protocol will retain
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

4.	In the md section, the technical options for *GROMACS* have to be specified:
    *	gmx_bin: The *GROMACS* command to perform mdrun
    *	mdp_names: The name of the mdp files present in the mdp folders specified above
    *	ngpus: The number of GPUs available
    *	mpi_procs: typically 1
    *	omp_procs: The number of threads used for each MD runs
    *	pinoffsets: Pinning the threads to specific positions to maximize the performance. This values depend on the GPU architecture
    *	use_tleap: True, this option is specified only if tleap is used to build the topology.

.. code-block:: console

    md:
        gmx_bin: gmx mdrun
        mdp_names:
            min_mdp: min.mdp
            nvt_mdp: nvt.mdp
            npt_mdp: npt.mdp
        ngpus: 4
        mpi_procs: 1
        omp_procs: 8
        pinoffsets: [ 0, 32, 64, 96 ]
        use_tleap: true

5.	In the target section, the chain ID of the target has to be specified.

.. code-block:: console

    target:
        chainID: [A]

7.	In the binder section, where the single point mutation will be performed, the following options have to be specified:
    *	chainID: The chain IDs of the binder
    *	mutating_chainID: The chain IDs where the mutation is desired. Since the mutating sequence is listed separately for each CDRs, the chain ID has to be a list also.
    *	mutating_resSeq: The residue sequence of the desired mutation sites. In this example, the sequences are listed separately for each CDRs.
    *	mutating_resname: The amino acid residues in one letter format, that correspond to the mutating_resSeq

.. code-block:: console

    binder:
        chainID: [B]
        mutating_chainID: [B,B,B]
        mutating_resSeq: [[220,221,222,223,224,225,226,227],[248,249,250,251,252,253,254],[294, 295, 296, 297, 298, 299, 300]]
        mutating_resname: [[S,G,F,D,F,S,D,A],[R,S,G,L,A,T,S],[K,S,R,R,G,Q,G]]


8.	In the scoring section, the choice of scoring functions have to be specified:
    *	functions: the choice of scoring functions, in this example, we use bluuesbmf, piepisa, evoef2, and MMPBSA
    *	consensus_threshold: The consensus threshold will be the criteria to decide whether the mutation is accepted or not. In this example, we set it to 3.
    *	nthreads: corresponds to the number of threads used to calculate scoring function
    *	mpiprocs: allows the MPI run for GMX_MMPBSA

.. code-block:: console

    scoring:
        functions: [ bluuesbmf, piepisa, evoef2, gmx_mmpbsa ]
        consensus_threshold: 3
        nthreads: 80
        mpiprocs: 2

Running the protocol
------------------------
Once the input file has been specified, and all the files are gathered, the protocol can now be run
by firstly activating the environment, if you haven't already.

.. code-block:: console

    mamba activate locuaz
    python /home/user/locuaz/locuaz/protocol.py config_tleap.yaml


Now the protocol will create the working directory folder. In this folder, the progress of the protocol is written
in the nb.log file. Afterwards, folders corresponding to each epochs and iterations will be created in this directory.

