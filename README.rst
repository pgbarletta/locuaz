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

Download mambaforge from:

https://github.com/conda-forge/miniforge

if on POWER9, go to the bottom, or just search 'POWER9'.

Can't install biobb using conda because It'll try to get packages that are not ppc64 compatible,
eg: gromacs, and some old babel version.

If on biobb 3.7, apply these patches:
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
    -On biobb_md/gromacs/pdb2gmx.py, line ~129:
        add:
            if self.dev:
                fu.log(f'Adding development options: {self.dev} -- DALE BOOO', self.out_log)
                self.cmd += self.dev.split()


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
