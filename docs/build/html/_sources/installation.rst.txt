.. highlight:: shell

============
Installation
============

Prerequisites
---------------

Given that pure conda is too slow, `Mambaforge <https://github.com/conda-forge/miniforge>`_ is
recommended instead.

Stable release
--------------

conda (recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: console

    mamba install locuaz

A drawback from the conda install is that when installing **biobb** through conda, it comes with some heavy
dependencies, like **GROMACS** itself.

pip
^^^

Create a conda environment from the :download:`usr_deps.yaml<../../usr_deps.yaml>`, which looks like this:

.. code-block:: console

    name: locuaz
    channels:
      - conda-forge
    dependencies:
      - conda-forge::python>=3.9,<3.10
      - conda-forge::ambertools>=22.0.0
      - conda-forge::tensorflow

Then, activate the environment and install the protocol through pip:

.. code-block:: console

    pip install locuaz

This is more involved that installing through conda, but the resulting environment won't be as heavy.

From sources
------------

Clone the `repo`_:

.. code-block:: console

    git clone https://github.com/pgbarletta/locuaz

And create the environment and install all the necessary dependencies at once:

.. code-block:: console

    mamba env create -f dev_deps.yaml

That's it. You can also change the environment's name by editing the `name` field of the `dev_deps.yml` file, before creating it.

Post-installation
------------------

If you want to use the ``dlp`` mutator You'll also have to get DLPacker's `weights <https://drive.google.com/file/d/1J4fV9aAr2nssrWN8mQ7Ui-9PVQseE0LQ/view?usp=sharing>`_
and place them on a dedicated ``dlpacker`` (actual name doesn't matter) directory, more info on :ref:`mutators:Mutators`.


.. _repo: https://github.com/pgbarletta/locuaz

For more info, check the :ref:`scoringfunctions:Scoring Functions`