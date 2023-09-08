===================================================
Tutorial: optimizing an antibody against a ligand
===================================================

Coming soon
=============



``allowed_nonstandard_residues``

remember that for scoring, the target chains are named A, and the ones from the binder, B.
write about the sanitization the protocol does with the splitted frames. All scorers but gmxmmpbsa use this PDBs.


``memory_positions`` and ``failed_memory_positions``:
empty memory slots on input user memory are allowed.
This allows the user to control for how many epochs will the non-empty memory be recalled.
Place them after the desired positions:

``memory_positions: [[2, 3, 4, 6, 7, 8], [], [], [] ]``


misc
""""

Add chainID with parmed::

        parm pre_nb.prmtop
        addPDB start.pdb
        outparm nb.prmtop
        go

You can also use cpptraj to add the chainIDs to the ``.rst7`` and ``.pdb`` files::

    parm nb.prmtop
    trajin nb.rst7
    trajout nb.pdb


asa

.. figure:: ./resources/tleap_iterations_dag.png
        :alt: iterations_dag

        Figure 1: snapshot of one optimized complex. **p53** is the yellow one on the left, with its loops colored red and
        blue, these loops have to be stabilized so it doesn't loose its function; the zinc atom and its coordinating
