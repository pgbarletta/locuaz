====================
Scorers
====================

All scorers (SFs) should be inside the ``['paths']['scorers']`` (see input config yaml) directory.
Their folder names should match the exact names of the scorers used in the config file and their binaries
should be on the top level of their folders and also be named with the exact scorer's name.
Some scorers have additional requirements, like parameter files, or the case of **gmxmmpbsa**
which is included with the protocol and only needs an input text file.

Additional requirements for specific SFs
---------------------------------------------


Assuming a scorers folder set to: ``['paths']['scorers']=home/user/my_SFs``.

gmxmmpbsa
^^^^^^^^^
| gmxmmpbsa directory: ``home/user/my_SFs/gmxmmpbsa``

This is the only scorer that comes bundled with the protocol.
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


We advise ``igb=2`` since ``igb=5`` requires a modification of the topology
that is not available in the protocol. From the Amber manual:

*With this option (igb=5), you should use the command "set default PBradii mbondi2"
in setting up the prmtop file, although "set default PBradii bondi" is also OK*

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

Though *locuaz*, as a rule, avoids stepping over other external tools, you'll notice that the input script
has 1-indexed ``startframe`` and ``endframe`` options that conflict with the 0-indexed ``start`` and ``end``
options from the input config file. This means that they have to be manually checked. If, for example, you are
writing 250 frames during your NPT run and you set your ``start`` option to ``50``, then your *gmxmmpbsa*
input file would look like the one above. Check :ref:`configurationfile:YAML configuration file` for more info.

pisa
^^^^^^
| pisa directory: ``home/user/my_SFs/pisa``
| pisa binary: ``home/user/my_SFs/pisa/pisa``
| pisa parameters: ``home/user/my_SFs/pisa/pisa.params``

rosetta
^^^^^^^^^
| rosetta directory: ``home/user/my_SFs/rosetta``

Symbolic links on the top rosetta folder should be added, pointing to files in the rosetta installation
Eg: inside the main rosetta folder, with the rosetta directory called **sources**::

    ln -s sources/rosetta_source/bin/InterfaceAnalyzer.linuxgccrelease rosetta
    ln -s sources/rosetta_database/ rosetta_database
    ln -s sources/rosetta_source/build/src/release/linux/4.14/64/ppc64le/gcc/8.4/ parameters
    ln -s sources/rosetta_source/build/external/release/linux/4.14/64/ppc64le/gcc/8.4/ external_parameters

haddock
^^^^^^^^^
| haddock directory: ``home/user/my_SFs/haddock``

As with all the scorers, all the necessary files have to be at the top level.
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
^^^^^^^^^
| piepisa directory: ``home/user/my_SFs/piepisa``

Download `pie <https://clsbweb.oden.utexas.edu/dock_details.html>`_. If you can run the binary, good,
if you can't, then you probably won't be able to run it, since compiling and running it in a
modern PC is quite cumbersome. Then, normalize the directory to the scorers standard:

* rename the **pie** folder to **piepisa**
* be sure to also have the **pisa** scorer
* Inside the **piepisa** folder, make symbolic links to the binaries and parameters so they have proper names::

    ln -s bin/pie_score pie
    ln -s bin/pie.params pie.params
    ln -s ../pisa/pisaEnergy_linux pisa
    ln -s ../pisa/pisa.params pisa.params

evoef2
^^^^^^^^^
| evoef2 directory: ``home/user/my_SFs/evoef2``

Download and compile `evoef2 <https://github.com/tommyhuangthu/EvoEF2>`_.

* rename the **EvoEF2** folder to **evoef2**
* Inside the **evoef2** folder, make a symbolic link to the binary so it has a proper name::

    ln -s bin/evoef2 evoef2

bluues
^^^^^^^
| bluues directory: ``home/user/my_SFs/bluues``

* Inside the **bluues** folder, make symbolic links to the binaries so it has a proper name::

    ln -s bin/bluues_new_2 bluues

bluuesbmf
^^^^^^^^^
| bluuesbmf directory: ``home/user/my_SFs/bluuesbmf``

* Inside the **bluuesbmf** folder, make symbolic links to the binary so it has a proper name::

    ln -s bin/bluues_new_2 bluues
    ln -s bin/score_bmf_3 bmf

autodockvina
^^^^^^^^^^^^^
| autodockvina directory: ``home/user/my_SFs/autodockvina``

Download `autodockvina <https://github.com/ccsb-scripps/AutoDock-Vina/releases>`_.
Then, normalize the directory to the scorers standard:
* create a folder named **autodockvina** with the downloaded binary
* Inside the **autodockvina** folder, make symbolic links to the binary so it has a proper name::

    ln -s vina_1.2.3_linux_x86_64 autodockvina



locuaz.scorers module
------------------------------

.. automodule:: locuaz.scorers
   :members:
   :undoc-members:
   :show-inheritance:

locuaz.abstractscorer module
-------------------------------------

.. automodule:: locuaz.abstractscorer
   :members:
   :undoc-members:
   :show-inheritance:

locuaz.bach
---------------------------------------

.. automodule:: locuaz.bach
   :members:
   :undoc-members:
   :show-inheritance:

locuaz.bluues
---------------------------------------

.. automodule:: locuaz.bluues
   :members:
   :undoc-members:
   :show-inheritance:

locuaz.bluuesbmf
---------------------------------------

.. automodule:: locuaz.bluuesbmf
   :members:
   :undoc-members:
   :show-inheritance:

locuaz.evoef2
---------------------------------------

.. automodule:: locuaz.evoef2
   :members:
   :undoc-members:
   :show-inheritance:

locuaz.haddock
---------------------------------------

.. automodule:: locuaz.haddock
   :members:
   :undoc-members:
   :show-inheritance:

locuaz.pisa
---------------------------------------

.. automodule:: locuaz.pisa
   :members:
   :undoc-members:
   :show-inheritance:

locuaz.piepisa
---------------------------------------

.. automodule:: locuaz.piepisa
   :members:
   :undoc-members:
   :show-inheritance:

locuaz.rosetta
---------------------------------------

.. automodule:: locuaz.rosetta
   :members:
   :undoc-members:
   :show-inheritance:

locuaz.gmxmmpbsa
---------------------------------------

.. automodule:: locuaz.gmxmmpbsa
   :members:
   :undoc-members:
   :show-inheritance:

locuaz.autodockvina
---------------------------------------

.. automodule:: locuaz.autodockvina
   :members:
   :undoc-members:
   :show-inheritance:

locuaz.scoring module
---------------------

.. automodule:: locuaz.scoring
   :members:
   :undoc-members:
   :show-inheritance:

locuaz.utils\_scoring module
----------------------------

.. automodule:: locuaz.utils_scoring
   :members:
   :undoc-members:
   :show-inheritance:
