Mutation Generators (deprecated)
================================

locuaz.abstractmutationgenerator
---------------------------------------

.. automodule:: locuaz.abstractmutationgenerator
   :members:
   :undoc-members:
   :show-inheritance:

locuaz.spm4
------------------

.. automodule:: locuaz.spm4
   :members:
   :undoc-members:
   :show-inheritance:

locuaz.spm4i
-------------------

.. automodule:: locuaz.spm4i
   :members:
   :undoc-members:
   :show-inheritance:

locuaz.spm4mmpbsa
------------------------

- ``gmxmmpbsa`` based generators like ``SPM4gmxmmpbsa`` need a residue decomposition file from ``gmxmmpbsa``,
  so the **gmxmmpbsa** script needs to include something along the lines of::

    /
    &decomp
    idecomp=2, dec_verbose=0,
    print_res="within 4"
    /

.. automodule:: locuaz.spm4mmpbsa
   :members:
   :undoc-members:
   :show-inheritance:
