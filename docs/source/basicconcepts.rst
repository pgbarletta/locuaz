==========================================
Basic concepts
==========================================

Main idea
-------------

locuaz is a protocol for the *in silico* optimization of antibodies and nanobodies.
The procedure begins with a random mutation in the binder sequence. Subsequently,
the bound conformations are sampled through molecular dynamics simulations, and the target and binder interactions
are assessed using various scoring functions. Finally, a consensus criterion is applied to the binding scores
in order to accept or reject the mutation. This process is repeated iteratively to explore new sequences
with potentially improved affinities towards their targets. If the mutation does not significantly improve the affinity,
then the mutants are discarded and a new set of mutants are generated, based on the original complex(es).
This workflow is outlined in Figure 1.

.. figure:: ./resources/protocol_workflow_simple.png
        :alt: workflow

        Figure 1: The protocol's workflow.

Each new set of mutants (or complexes) generated will be referred to as an **epoch**, while each complex is referred to
as an **iteration**. 

Units
------------

locuaz has to coordinate between several external programs and be flexible enough to allow different
protocols to be run, hence, some abstractions are needed. We will call these abstractions *units*.

Throughout this document, we will refer to the user configuration options as ``config``, and its various options as
``config["main"]["name"]``, ``config["scoring"]["functions"]``, etc...

Mutation Generator
^^^^^^^^^^^^^^^^^^^^
These units are the one in charge of generating the new binders. These are the currently available generators:

SPM4
"""""
This is a Single Point Mutation generator. This means that it chooses a single position (from the user input
``config["binder"]["mutating_resSeq"]``), and all the mutations will be performed there.
To choose which amino acid will be used, it splits all amino acids (except cysteine, which is discarded) in the
following categories: **negative**, **positive**, **hydrophobic** and **ring-containing**.
Then, it chooses 1 from each group to generate as many mutations as the user asked for
(``config["protocol"]["branches"]`` option).

Set ``config["generation"]["generator"]`` to ``SPM4`` use this generator.

SPM4i
"""""
Same as ``SPM4``, but not any position in ``config["binder"]["mutating_resSeq"]`` may be mutated. Those that are not
part of the current interface will be discarded. To determine the interface, locuaz uses the **freesasa** library which
uses a rolling-probe, whose radius can be set using the ``config["generation"]["probe_radius"]`` to any value ranging
from ``0.1`` to ``4.0`` (in angstrom units). The bigger the radius, the more residues will be classified as part of
the interface; the default is ``1.4``.

Set ``config["generation"]["generator"]`` to ``SPM4i`` use this generator.

SPM4gmxmmpbsa
"""""""""""""""
Same as ``SPM4i``, but besides **freesasa**, it's based on the use as **gmxmmpbsa** scoring function. The generator
will read the **decomp_gmxmmpbsa.csv** output file from **gmxmmpbsa** and pick the residue that is collaborating the
list with the interaction with the target. Obviously, this position has to also comply with the previous prerequisites,
that is, being part of the interface and one of the positions included in  ``config["binder"]["mutating_resSeq"]``.
You can also set the probe radius in this mutator.

Don't forget to include ``gmxmmpbsa`` alongside your other scoring functions (in ``config["scoring"]["functions"]``),
and to include instructions in the **gmxmmpbsa** input file to perform the decompositions. The decomposition section
should look something like this:

.. code-block:: console

    /
    &decomp
    idecomp=2, dec_verbose=0,
    print_res="within 4"
    /

Check Amber's manual and `gmx-mmpbsa <https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/input_file/>`_ docs for more info.

Set ``config["generation"]["generator"]`` to ``SPM4gmxmmpbsa`` use this generator.

Mutator
^^^^^^^^^^
T

evoef2
"""""""
ss

MD
^^^^^^^^^^
as


Summary
--------

The

.. figure:: ./resources/protocol_workflow.png
        :alt: enhanced workflow

        Figure 2: the protocol's main concepts and the stages at which they act. An **iteration** is highlighted in green
        and the **epoch** in pink.

dfd-