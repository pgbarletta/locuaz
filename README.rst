====
locuaz
====


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

* Not tested for multiple chains


Install
------

Download Mambaforge from:

https://github.com/conda-forge/miniforge

if on POWER9, go to the bottom, or just search 'POWER9'.

Can't install biobb using conda because It'll try to get packages that are not ppc64 compatible,
eg: gromacs, and some old babel version.

* Make sure your cmake version is >3.17. If on MARCONI100:

```
module load autoload cmake
```

If on MDAnalysis 2.24.2 or older:
- On MDAnalysis/topology/tpr/utils.py line 330:
```
segid = f"seg_{i}_{molblock}"
```

    replace with:

```
segid = molblock[14:] if molblock[:14] == "Protein_chain_" else molblock
```



On scoring
--------

All scoring functions (SFs) should be inside the config['paths']['scoring_functions'] directory.
Their folder names should match the exact SF names used in the config file and their binaries
should be on the top level of their folders and also be named with the exact SF name. 
Eg: for the pisa scoring function and config['paths']['scoring_functions']='home/user/my_SFs'

pisa directory: 'home/user/my_SFs/pisa'
pisa binary: 'home/user/my_SFs/pisa/pisa'
pisa parameters (pisa's a special case): 'home/user/my_SFs/pisa/pisa.params'

Additional requirements for specific SFs:
 - pisa: see above.
 - rosetta: symbolic links on the top rosetta folder should be added, pointing the InterfaceAnalyzer,
   the database, the parameters directory and the external parameters directory. 
   Eg: inside the main rosetta folder
    ln -s sources/rosetta_source/bin/InterfaceAnalyzer.linuxgccrelease rosetta
    ln -s sources/rosetta_database/ rosetta_database
    ln -s sources/rosetta_source/build/src/release/linux/4.14/64/ppc64le/gcc/8.4/ parameters
    ln -s sources/rosetta_source/build/external/release/linux/4.14/64/ppc64le/gcc/8.4 external_parameters

 - haddock:
    The 'template_scoring.inp' has to be at the top level of the haddock
    The 'rescoring-scripts' folder has to be at the top level of the haddock
    The 'haddock' folder has to be at the top level of the haddock
    The 'cns_solve' or 'cns_solve_X.Y' (where 'X'.'Y' is the version number) folder
        has to be at the top level of the haddock

    ln -s ./cns_solve_1.3/ibm-ppc64le-linux/bin/cns cns
    ln -s haddock/protocols protocols
    ln -s haddock/toppar/ toppar
    ln -s cns_solve_1.3/cns_solve_env cns_solve_env
    ln -s haddock/haddock_configure.csh haddock_configure.csh

Features
--------

 - If you want to use amber topologies:

```
mamba install ambertools 
```


Credits
-------

- Biobb:
    https://mmb.irbbarcelona.org/biobb/documentation/source
    https://mmb.irbbarcelona.org/biobb/workflows/tutorials/md_setup


