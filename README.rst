========
locuaz
========


.. image:: https://img.shields.io/pypi/v/locuaz.svg
        :target: https://pypi.python.org/pypi/locuaz

.. image:: https://img.shields.io/travis/pgbarletta/locuaz.svg
        :target: https://app.travis-ci.com/github/pgbarletta/locuaz/builds

.. image:: https://readthedocs.org/projects/locuaz/badge/?version=latest
        :target: https://locuaz.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status


Looping Uniquely Catered Amino Acid Sequences


* Free software: MIT license
* Documentation: https://locuaz.readthedocs.io.

Install
--------

Mambaforge is recommended instead of pure conda. Download Mambaforge from:

https://github.com/conda-forge/miniforge

Clone this repo and, optionally, get the **DLPacker**  submodule as well::

    git clone https://github.com/pgbarletta/locuaz
    git submodule int
    git submodule update

You'll also have to get DLPacker's `weights <https://drive.google.com/file/d/1J4fV9aAr2nssrWN8mQ7Ui-9PVQseE0LQ/view?usp=sharing>`_.

Post-Install
-------------
If on MDAnalysis 2.4.3 or older, edit the file ``MDAnalysis/topology/tpr/utils.py`` line 330::
    
  segid = f"seg_{i}_{molblock}"

replace it with::

    segid = molblock[14:] if molblock[:14] == "Protein_chain_" else molblock


On scoring
----------------

All scoring functions (SFs) should be inside the ``['paths']['scoring_functions']`` (see input config yaml) directory.
Their folder names should match the exact SF names used in the config file and their binaries
should be on the top level of their folders and also be named with the exact SF name.
Some scoring functions have additional requirements, like parameter files,
or the case of **gmxmmpbsa** which is included with the protocol and only needs an input text file.

Additional requirements for specific SFs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Assuming a scoring functions folder set to: ``['paths']['scoring_functions']=home/user/my_SFs``.

gmxmmpbsa
""""""""""
| gmxmmpbsa directory: ``home/user/my_SFs/gmxmmpbsa``

This is the only scoring function that comes bundled with the protocol.
Inside the **gmxmmpbsa** folder, a **gmxmmpbsa** input text file is needed.
The contents are up to the user. For example, for a simple MM-GBSA::

    &general
    sys_name="Prot-Prot",
    startframe=51,
    endframe=250,
    /
    &gb
    igb=2, saltcon=0.150,
    /

And if residue decomposition is needed (for a mmpbsa generator)::

    &general
    sys_name="Prot-Prot",
    startframe=51,
    endframe=250,
    /
    &gb
    igb=2, saltcon=0.150,
    /
    &decomp
    idecomp=2, dec_verbose=0,
    print_res="within 4"
    /

pisa
"""""
| pisa directory: ``home/user/my_SFs/pisa``
| pisa binary: ``home/user/my_SFs/pisa/pisa``
| pisa parameters: ``home/user/my_SFs/pisa/pisa.params``

rosetta
"""""""""
| rosetta directory: ``home/user/my_SFs/rosetta``

Symbolic links on the top rosetta folder should be added, pointing to files in the rosetta installation
Eg: inside the main rosetta folder, with the rosetta directory called **sources**::

    ln -s sources/rosetta_source/bin/InterfaceAnalyzer.linuxgccrelease rosetta
    ln -s sources/rosetta_database/ rosetta_database
    ln -s sources/rosetta_source/build/src/release/linux/4.14/64/ppc64le/gcc/8.4/ parameters
    ln -s sources/rosetta_source/build/external/release/linux/4.14/64/ppc64le/gcc/8.4/ external_parameters

haddock
""""""""
| haddock directory: ``home/user/my_SFs/haddock``

As with all the scoring functions, all the necessary files have to be at the top level.
The **template_scoring.inp** file has to be at the top level of the haddock, as the **rescoring-scripts** folder
(included with the protocol insed the **sample_bin** folder).
Then, the following smybolic links have to be created.
Version number and specific folder names and locations may change::

    ln -s ./cns_solve_1.3/ibm-ppc64le-linux/bin/cns haddock
    ln -s haddock/protocols/ protocols
    ln -s haddock/toppar/ toppar
    ln -s cns_solve_1.3/cns_solve_env cns_solve_env
    ln -s haddock/haddock_configure.csh haddock_configure.csh

piepisa
""""""""
| piepisa directory: ``home/user/my_SFs/piepisa``

Download `pie <https://clsbweb.oden.utexas.edu/dock_details.html>`_. If you can run the binary, good,
if you can't, then you probably won't be able to run it, since compiling and running it in a
modern PC is quite cumbersome. Then, normalize the directory to the scoring functions standard:

* rename the **pie** folder to **piepisa**
* be sure to also have the **pisa** scoring function
* Inside the **piepisa** folder, make symbolic links to the binaries and parameters so they have proper names::

    ln -s bin/pie_score pie
    ln -s bin/pie.params pie.params
    ln -s ../pisa/pisaEnergy_linux pisa
    ln -s ../pisa/pisa.params pisa.params

evoef2
""""""
| evoef2 directory: ``home/user/my_SFs/evoef2``

Download and compile `evoef2 <https://github.com/tommyhuangthu/EvoEF2>`_.

* rename the **EvoEF2** folder to **evoef2**
* Inside the **evoef2** folder, make a symbolic link to the binary so it has a proper name::

    ln -s bin/evoef2 evoef2

bluues
""""""""
| bluues directory: ``home/user/my_SFs/bluues``

* Inside the **bluues** folder, make symbolic links to the binaries so it has a proper name::

    ln -s bin/bluues_new_2 bluues

bluuesbmf
"""""""""
| bluuesbmf directory: ``home/user/my_SFs/bluuesbmf``

* Inside the **bluuesbmf** folder, make symbolic links to the binary so it has a proper name::

    ln -s bin/bluues_new_2 bluues
    ln -s bin/score_bmf_3 bmf

autodockvina
""""""""""""
| autodockvina directory: ``home/user/my_SFs/autodockvina``

Download `autodockvina <https://github.com/ccsb-scripps/AutoDock-Vina/releases>`_.
Then, normalize the directory to the scoring functions standard:
* create a folder named **autodockvina** with the downloaded binary
* Inside the **autodockvina** folder, make symbolic links to the binary so it has a proper name::

    ln -s vina_1.2.3_linux_x86_64 autodockvina

Mutators
---------

-  DLPacker is included as a submodule. To download it::

    git submodule init
    git submodule update

Then, in a ``dlpacker`` directory, the following files have to be present:

1. ``charges.rtp``
2. ``library.npz``
3. ``DLPacker_weights.h5``

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
