.. highlight:: shell

============
Installation
============

There are 3 options for installing locuaz. Conda + pip (Option 1) is recommended, since its both lightweight
and quick. A pure conda install is also available (Option2). It is also possible to compile from source (Option 3).
In all cases a :ref:`installation:Post-installation` step is required.

Prerequisites
---------------

We recommend using `Miniforge <https://github.com/conda-forge/miniforge>`_ instead of pure conda,
since it's much faster. We used to recommend miniforge but it's been
`discontinued <https://conda-forge.org/news/2024/07/29/sunsetting-mambaforge/>`_.

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
      - conda-forge::python>=3.10
      - conda-forge::ambertools>=22.0.0
      - conda-forge::tensorflow
      - conda-forge::openbabel

To create the environment:

.. code-block:: console

    mamba env create -f usr_deps.yaml


Then, activate the environment and install the protocol through pip:

.. code-block:: console

    mamba activate locuaz
    pip install locuaz

This option takes an extra step with respect to using only conda, but the process will be faster and the
resulting environment won't be as heavy.

apptainer container (Option 2)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
You can pull any version from the GitHub repository like this:

.. code-block:: console

    apptainer pull oras://ghcr.io/pgbarletta/locuaz.sif:0.6.2

These containers will be signed with a public key, so you can verify them::

    apptainer verify --url https://keys.openpgp.org locuaz.sif

*apptainer* containers are available for a plug-and-play approach. Just change
the way you call *locuaz*. If you'd call the *locuaz* entry-point like this:

.. code-block:: console

    locuaz config.yaml

this is how you'd call it when inside a container:

.. code-block:: console

    apptainer exec --nv locuaz.sif locuaz config.yaml

Check :ref:`jobsubmission:Job submission` for more info on using *locuaz* with
containers.

From sources (Option 3, developers only)
-----------------------------------------

Clone the `repo`_:

.. code-block:: console

    git clone https://github.com/pgbarletta/locuaz

Create the environment and install all the necessary dependencies at once:

.. code-block:: console

    mamba env create -f dev_deps.yaml

And inside the newly cloned dir, install *locuaz* in development mode: ::

    pip install -e .

That's it. You can also change the environment's name by editing the `name` field of the `dev_deps.yml` file, before creating it.

Post-installation
------------------

If you want to use the ``dlp`` mutator You'll also have to get DLPacker's `weights`_ and place them on a dedicated
``dlpacker`` (actual name doesn't matter) directory, more info on :ref:`mutators:Mutators`.
You also may want extra scorers, check the :ref:`scorers:scorers`

.. _repo: https://github.com/pgbarletta/locuaz
.. _weights: https://drive.google.com/file/d/1J4fV9aAr2nssrWN8mQ7Ui-9PVQseE0LQ/view?usp=sharing
