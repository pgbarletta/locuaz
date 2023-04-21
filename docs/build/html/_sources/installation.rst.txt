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

conda (mamba)
^^^^^^^^^^^^^^

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
      - conda-forge::python>=3.9.15,<=3.10.10
      - conda-forge::ambertools>=22.0.0
      - conda-forge::tensorflow

Then, activate the environment and install the protocol through pip:

.. code-block:: console

    pip install locuaz

This is more involved that installing through conda, but the resulting environment won't be as heavy.

Why there's no straight pip install
""""""""""""""""""""""""""""""""""""
Being this a high-level protocol it has many dependencies and, unfortunately, python packaging has its quirks, so
different developers have solved their issues in different ways.
A straight `pip install` is out of the question given the tensorflow (for DLPacker) and ambertools (for tleap)
dependencies which are only easily available through conda, due to their involved installation process.


From sources
------------

Clone the `repo`_ and, optionally, get the **DLPacker**  submodule as well:

.. code-block:: console

    git clone https://github.com/pgbarletta/locuaz
    git submodule int
    git submodule update

Finally, create the environment and install all the necessary dependencies at once:

    mamba env create -f dev_deps.yaml

That's it. You can also change the environment's name by editing the `name` field of the `dev_deps.yml` file, before creating it.

Post-installation
------------------

The protocol
You'll also have to get DLPacker's `weights <https://drive.google.com/file/d/1J4fV9aAr2nssrWN8mQ7Ui-9PVQseE0LQ/view?usp=sharing>`_
and place them on a dedicated ``dlpacker`` (actual name doesn't matter) directory, more info on the dedicated Mutators section.


.. _repo: https://github.com/pgbarletta/locuaz

For more info, check the :ref:`Scoring Functions`