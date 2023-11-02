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

Create a conda environment YAML file named, for example, ``usr_deps.yaml``::

    name: locuaz
    channels:
      - conda-forge
    dependencies:
      - conda-forge::python>=3.10
      - conda-forge::ambertools>=22.0.0
      - conda-forge::tensorflow
      - conda-forge::openbabel
      - conda-forge::pygraphviz

by running:

.. code-block:: console

    mamba env create -f usr_deps.yaml

Then, activate the environment and install the protocol through pip:

.. code-block:: console

    pip install locuaz

Check the `Installation`_ section on the docs for more info.


Post-Install
-------------

Using DLPacker based mutators
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you wish to used *DLPacker* based *Mutators* (``dlp`` and ``dlpr``), DLPacker weights have to be downloaded.
Get them `here <https://drive.google.com/file/d/1J4fV9aAr2nssrWN8mQ7Ui-9PVQseE0LQ/view?usp=sharing>`_
or `here`_. These weights have to be extracted to a dedicated folder and its path has to be specified in the
input config under the ``paths`` key, on the  ``mutator`` option. Check the `docs`_ for more info.

Citing
-------



Credits
-------

- `Biobb <https://mmb.irbbarcelona.org/biobb/documentation/source>`_
- `MDAnalysis <https://github.com/MDAnalysis/mdanalysis>`_
- `FreeSASA <https://github.com/freesasa/freesasa-python>`_

.. _docs: https://locuaz.readthedocs.io/en/latest/mutators.html
.. _Installation: https://locuaz.readthedocs.io/en/latest/installation.html
