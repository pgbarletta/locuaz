LOCUAZ optimization protocol
======================================

.. image:: https://img.shields.io/pypi/v/locuaz.svg
        :target: https://pypi.python.org/pypi/locuaz

.. image:: https://readthedocs.org/projects/locuaz/badge/?version=latest
        :target: https://locuaz.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status

.. image:: https://img.shields.io/badge/cite-locuaz-red
        :target: citing.html
        :alt: Cite LOCUAZ

.. image:: https://img.shields.io/badge/license-MIT-yellow
        :target: citing.html
        :alt: License

*locuaz* is a protocol for the *in silico* optimization of antibodies and antibody fragments, such as nanobodies.
It could also be employed to optimise other objects such as (poli)peptides or other binders.

The procedure begins with random mutations in the binder sequence, which generate different target-binder
complexes. Subsequently, the complexes are minimized, equilibrated with a NVT run and then sampled
through a NPT Molecular Dynamics (MD) simulation, so the target and binder interactions can be assessed
with the chosen scoring method(s)
Finally, a consensus criterion is applied to the binding scores in order to accept or reject the mutation.
If the mutation does not significantly improve the affinity, then the mutants are discarded and a new
set of mutants are generated, based on the original complex(es). On the other hand, when accepted, the
complex is fed into the next run of the process, to keep exploring new sequences with potentially improved affinities
towards their targets.

This workflow is outlined in Figure 1.

.. figure:: ./resources/protocol_workflow_simple.png
        :alt: workflow
        :scale: 75%

        Figure 1: The protocol's workflow.



.. toctree::
    :maxdepth: 2
    :caption: Installation

    installation

.. toctree::
    :maxdepth: 2
    :caption: Learning

    basicconcepts
    tutorialsimple
    tutorialtleap
    tutorialligand
    jobsubmission

.. toctree::
    :maxdepth: 1
    :caption: Reference

    mutationgenerators
    scoringfunctions
    mutators
    statistics
    pruners
    misc
    configurationfile
    citing
    history

Indices and tables
==================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
