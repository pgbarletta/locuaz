===================================================
Frequently Asked Questions
===================================================

.. _faq1:

1\_ |q1|
--------
*tleap* is a bit of an old piece of software and it seems that back then chainID
info wasn't very important. If you write out a PDB file along with your *prmtop*
and *rst7* files when building your topology, you'll see that the PDB has no
chainID info and neither does the topology.

Fortunately, *parmed* ---which will be included with your locuaz installation---,
added some new flags to the ``.prmtop`` format, including ``RESIDUE_CHAINID``,
which will hold the information we want. In order to add it, along with some other
extra info, we will ask parmed to get it from our initial PDB, which should have
chainID info

This short parmed script will do it::

        parm no_chainid_top.prmtop
        addPDB start.pdb
        outparm top.prmtop
        go

``no_chainid_top.prmtop`` is the topology we got out from tleap, with water, ions
and all that. ``start.pdb`` is the PDB tleap started out from, without water, box
or anything, but with the chainID information.

We then run this script::

    parmed -i add_chainID.parmed

You can now check ``top.prmtop`` for parmed additions::

    %FLAG RESIDUE_CHAINID
    %COMMENT Residue chain ID (chainId) read from PDB file; DIMENSION(NRES)
    %FORMAT(20a4)
    A   A   A   A   A   A   A   A   A   A   A   A   A   A   B   B   B   B   B   B

The final step is to get a PDB out of this topology and the *rst7* tleap gave us.
We can, again, write a parmed script for this::

    parm top.prmtop
    trajin rst7_from_tleap.rst7
    trajout pdb_with_chainID.pdb

Or, we can use another program that comes with *ambertools*, which was made for
this specific purpose::

    ambpdb -p top.prmtop -c rst7_from_tleap.rst7 > pdb_with_chainID.pdb

Now we can use ``pdb_with_chainID.pdb`` to start out our optimization, probably
after giving it a better name.

.. |q1| replace:: I used tleap to build my initial system
    but now my starting PDB has no chainID information and locuaz won't take it.
    What can I do about it?


.. _faq2:

2\_ My complex is spreading apart. Can I use GROMACS pulling to keep them together?
-----------------------------------------------------------------------------------
Of course! This GROMACS feature is orthogonal to the optimization so *locuaz*
doesn't concern itself with it.
What it will do is give GROMACS an index selection
file (``.ndx``) with selections for the target and the binder. You can then
reference **target** and **binder** as the pull groups in the pull code.
This pull code goes in the mdp file you'll use for the NPT run::


    ; Pull code
    pull                        = yes
    pull-pbc-ref-prev-step-com  = yes       ;
    pull-ncoords                = 1         ; only one reaction coordinate
    pull-ngroups                = 2         ; two groups defining one reaction coordinate
    pull-group1-name            = target
    pull-group1-pbcatom         = 4075      ; atom close to the center
    pull-group2-name            = binder
    pull-group2-pbcatom         = 7251      ; atom close to the center
    pull-coord1-groups          = 1 2       ; target and binder define the reaction coordinate
    pull-coord1-type            = umbrella  ; harmonic potential
    pull-coord1-geometry        = distance  ;
    pull-coord1-dim             = Y Y Y     ; pull along all cordinates
    pull-coord1-start           = yes       ; define initial COM distance > 0
    pull-coord1-rate            = 0.00      ; keep it fixed
    pull-coord1-k               = 1         ; kJ mol^-1 nm^-2


Notice the **target** and **binder** as the pull groups 1 and 2, respectively.
You can probably use this same sample file and just modify the atoms at
``pull-group1-pbcatom`` and ``pull-group2-pbcatom``, and if necessary, the
``pull-coord1-k``.

But do remember to lift this pulling and optimize your binder using free MD
once you have a somewhat better binder than can stay with its target without any
restraints. Try optimizing for 5 epochs and then lift the restrictions to keep
optimizing without pulling.

.. _faq3:

3\_ What is NUMA? How can I know the number of NUMA regions my machine has?

NUMA (Non-Unified Memory Access) it's a design idea to deal with multiple CPU
processors and their access to memory, and GPUs as well.

For us users, it's important to be NUMA-aware, that is, understand that some
cores have faster access to certain regions and slower access to other ones.
The same applies to the relationship between CPUs and GPUs.

To know your NUMA layout, just run this command::

    $ numactl  --hardware
    available: 2 nodes (0-1)
    node 0 cpus: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
    node 0 size: 256924 MB
    node 0 free: 249438 MB
    node 1 cpus: 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
    node 1 size: 257999 MB
    node 1 free: 256694 MB
    node distances:
    node   0   1
      0:  10  11
      1:  11  10

On this cluster node, there are 2 NUMA regions (nodes), with 16 cores (32 threads)
each. I'll add that these nodes correspond to each of the physical CPU sockets.
By the distances, you can see that it's more convenient for the threads of each
socket to access the memory that belongs to that same node.

But for our purpose, GPU-accelerated MD, we don't care much about RAM memory,
but about GPU-CPU affinity, that is, to which threads the GPUs are closer.
Let's look at the next command to answer this::

    $ nvidia-smi topo -m
            GPU0    GPU1    GPU2    GPU3    mlx5_0  mlx5_1  mlx5_2  mlx5_3  CPU Affinity    NUMA Affinity
    GPU0     X      NV4     NV4     NV4     PXB     NODE    NODE    NODE    0-15    0
    GPU1    NV4      X      NV4     NV4     NODE    PXB     NODE    NODE    0-15    0
    GPU2    NV4     NV4      X      NV4     NODE    NODE    PXB     NODE    0-15    0
    GPU3    NV4     NV4     NV4      X      NODE    NODE    NODE    PXB     0-15    0
    mlx5_0  PXB     NODE    NODE    NODE     X      NODE    NODE    NODE
    mlx5_1  NODE    PXB     NODE    NODE    NODE     X      NODE    NODE
    mlx5_2  NODE    NODE    PXB     NODE    NODE    NODE     X      NODE
    mlx5_3  NODE    NODE    NODE    PXB     NODE    NODE    NODE     X

    Legend:

      X    = Self
      SYS  = Connection traversing PCIe as well as the SMP interconnect between NUMA nodes (e.g., QPI/UPI)
      NODE = Connection traversing PCIe as well as the interconnect between PCIe Host Bridges within a NUMA node
      PHB  = Connection traversing PCIe as well as a PCIe Host Bridge (typically the CPU)
      PXB  = Connection traversing multiple PCIe bridges (without traversing the PCIe Host Bridge)
      PIX  = Connection traversing at most a single PCIe bridge
      NV#  = Connection traversing a bonded set of # NVLinks

There is a lot of information on the output, but we only care about the ``CPU Affinity`` column.
This tells us that all 4 GPUs have a direct link to threads 0-15, which means
threads 16-31 will take longer to communicate with the GPUs.

This is actually quite unexpected. I'd expect GPUs 0 and 1 to be closer to threads
0-15 and GPUs 2 and 3 to be closer to threads 16-31. I'll contact the sysadmins
and hopefully I'll remember to update this site.


.. _faq4:

4\_ Why use both ``memory_positions`` and ``failed_memory_positions`` at the same time?

.. _faq5:

5\_ What are the empty brackets in the ``memory_positions`` list?

empty memory slots on input user memory are allowed.
This allows the user to control for how many epochs will the non-empty memory be recalled.
Place them after the desired positions:

``memory_positions: [[2, 3, 4, 6, 7, 8], [], [], [] ]``

