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

Download Mambaforge from:

https://github.com/conda-forge/miniforge

if on POWER9, go to the bottom, or just search 'POWER9'.

Can't install biobb using conda because It'll try to get packages that are not ppc64 compatible,
eg: gromacs, and some old babel version.

* Make sure your cmake version is >3.17. If on MARCONI:

```
module load autoload cmake
```

### If on biobb 3.8, apply these patches:
    - On biobb_analysis/gromacs/gmximage.py, line 126:
      replace:
```
# If fitting provided, echo fit_selection
if self.fit == 'none':
    if self.center:
        selections = '\"' + self.center_selection + '\" \"' + self.output_selection + '\"'
    else:
        selections = '\"' + self.output_selection + '\"'
else:
    if self.center:
        selections = '\"' + self.fit_selection + '\" \"' + self.center_selection + '\" \"' + self.output_selection + '\"'
    else:
        selections = '\"' + self.fit_selection + '\" \"' + self.output_selection + '\"'
```
        with:
```
# If fitting provided, echo fit_selection
if self.fit == 'none':
    selections = '\"' + self.center_selection + '\" \"' + self.output_selection + '\"'
else:
    selections = '\"' + self.fit_selection + '\" \"' + self.center_selection + '\" \"' + self.output_selection + '\"'
```

### If on biobb 3.7, apply these patches:
    - On biobb_md/gromacs/mdrun.py, line 201:
        replace:
            self.cmd += [self.dev.split()]
        with:
            self.cmd.append(self.dev)
    -On biobb_md/gromacs/solvate.py, line ~82:
        add:
            self.dev = properties.get('dev')
    -On biobb_md/gromacs/solvate.py, line ~129:
        add:
            if self.dev:
                fu.log(f'Adding development options: {self.dev} -- DALE BOOO', self.out_log)
                self.cmd += self.dev.split()
    -On biobb_md/gromacs/pdb2gmx.py, line ~82:
        add:
            self.dev = properties.get('dev')
    -On biobb_md/gromacs/pdb2gmx.py, line ~127:
        add:
            if self.dev:
                fu.log(f'Adding development options: {self.dev} -- DALE BOOO', self.out_log)
                self.cmd += self.dev.split()
    -On biobb_analysis/gromacs/gmx_trjconv_str_ens.py, line 82:
        replace:
            self.fit_selection = properties.get('fit_selection', "System")
        with:
            self.selection = properties.get('selection', "System")
    -On biobb_md/gromacs/editconf.py, line ~74:
        add:
            self.dev = properties.get('dev')
    -On biobb_md/gromacs/solvate.py, line ~117:
        add:
            if self.dev:
                fu.log(f'Adding development options: {self.dev} -- DALE BOOO', self.out_log)
                self.cmd += self.dev.split()


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

* TODO

Credits
-------

- Biobb:
    https://mmb.irbbarcelona.org/biobb/documentation/source
    https://mmb.irbbarcelona.org/biobb/workflows/tutorials/md_setup

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
