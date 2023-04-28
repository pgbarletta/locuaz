========
locuaz
========


.. image:: https://img.shields.io/pypi/v/locuaz.svg
        :target: https://pypi.python.org/pypi/locuaz

.. image:: https://readthedocs.org/projects/locuaz/badge/?version=latest
        :target: https://locuaz.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status


Looping Uniquely Catered Amino Acid Sequences


* Free software: MIT license
* Documentation: https://locuaz.readthedocs.io.

Install
--------

Check the `Installation <https://locuaz.readthedocs.io/en/latest/installation.html>`_ section on the docs.


Post-Install
-------------
If on MDAnalysis 2.4.3 or older, edit the file ``MDAnalysis/topology/tpr/utils.py`` line 330::
    
  segid = f"seg_{i}_{molblock}"

replace it with::

    segid = molblock[14:] if molblock[:14] == "Protein_chain_" else molblock


On scoring
----------------


Mutators
---------

-  DLPacker is included as a submodule. To download it::

    git submodule init
    git submodule update

Then, in

The first 2 can be copied from the recently downloaded directory (``locuaz/DLPacker``).
The weights have to be `downloaded <https://drive.google.com/file/d/1J4fV9aAr2nssrWN8mQ7Ui-9PVQseE0LQ/view?usp=sharing>`_.
Then, the path to the ``dlpacker`` directory has to be specified in the input config under the
``paths`` key, on the  ``mutator`` option.

Generators
-----------

- ``gmxmmpbsa`` based generators like ``SPM4gmxmmpbsa`` need a residue decomposition file from ``gmxmmpbsa``,
  so the **gmxmmpbsa** script needs to include something along the lines of::

    /
    &decomp
    idecomp=2, dec_verbose=0,
    print_res="within 4"
    /

Credits
-------

- `Biobb <https://mmb.irbbarcelona.org/biobb/documentation/source>`_
- `MDAnalysis <https://github.com/MDAnalysis/mdanalysis>`_
- `FreeSASA <https://github.com/freesasa/freesasa-python>`_
