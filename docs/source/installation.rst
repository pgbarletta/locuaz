.. highlight:: shell

============
Installation
============

There are 3 options for installing locuaz. Conda + pip (Option 1) is recommended, since its both lightweight
and quick. A pure conda install is also available (Option2). It is also possible to compile from source (Option 3).
In all cases a :ref:`installation:Post-installation` step is required.

Prerequisites
---------------

We recommend using `Mambaforge <https://github.com/conda-forge/miniforge>`_ instead of pure conda,
since it's much faster

Stable release
--------------

conda + pip (Option 1)
^^^^^^^^^^^^^^^^^^^^^^^^^

Create a conda environment from the :download:`usr_deps.yaml<../../usr_deps.yaml>`, which looks like this:

.. code-block:: console

    name: locuaz
    channels:
      - conda-forge
    dependencies:
      - conda-forge::python>=3.9,<3.10
      - conda-forge::ambertools>=22.0.0
      - conda-forge::tensorflow

To create the environment:

.. code-block:: console

    mamba create -f usr_deps.yaml


Then, activate the environment and install the protocol through pip:

.. code-block:: console

    mamba activate locuaz
    pip install locuaz

This option takes an extra step with respect to using only conda, but the process will be faster and the
resulting environment won't be as heavy.

conda (Option 2)
^^^^^^^^^^^^^^^^

.. code-block:: console

    mamba install locuaz

A drawback from the conda install is that when installing **biobb** through conda, it comes with some heavy
dependencies, like **GROMACS** itself.

From sources (Option 3)
------------------------

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