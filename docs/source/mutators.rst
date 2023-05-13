Mutators 
===========

Using DLPacker based mutators
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you wish to used *DLPacker* based *Mutators* (``dlp`` and ``dlpr``), DLPacker weights have to be downloaded.
Get them `over here <https://drive.google.com/file/d/1J4fV9aAr2nssrWN8mQ7Ui-9PVQseE0LQ/view?usp=sharing>`_
or ` here`_. These weights have to be extracted to a dedicated folder and its path has to be specified in the
input config under the ``paths`` key, on the  ``mutator`` option.

locuaz.mutation module
----------------------

.. automodule:: locuaz.mutation
   :members:
   :undoc-members:
   :show-inheritance:

locuaz.mutators module
----------------------

.. automodule:: locuaz.mutators
   :members:
   :undoc-members:
   :show-inheritance:

locuaz.basemutator module
-------------------------

.. automodule:: locuaz.basemutator
   :members:
   :undoc-members:
   :show-inheritance:

locuaz.mutatordlp module
------------------------

.. automodule:: locuaz.mutatordlp
   :members:
   :undoc-members:
   :show-inheritance:

locuaz.mutatordlpr module
-------------------------

The ``dlpr`` mutator uses *DLPacker* to reorient the side-chains of the residues that are within a certain distance
of the mutated position. To adjust this distance use the ``config['mutation']['reconstruct_radius']`` option.
Notice that non-standard residues and cysteins that are forming sulfide bridges will be excluded from this
reorientation.


.. automodule:: locuaz.mutatordlpr
   :members:
   :undoc-members:
   :show-inheritance:

locuaz.mutatorevoef2 module
---------------------------

.. automodule:: locuaz.mutatorevoef2
   :members:
   :undoc-members:
   :show-inheritance:


.. _here: https://istitutoitalianotecnologia-my.sharepoint.com/:u:/g/personal/walter_rocchia_iit_it/Efzdf2sgKwJNmJskcHDE7yUBQMVgFsbpACeQLDGRYKvQOA?e=2E0daX