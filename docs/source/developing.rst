===========
Developing
===========

Releasing a new version
------------------------
The following steps have to be completed in order to release a new version of *locuaz*:

1. Make sure you have all the dependencies::

    pip install twine build sphinx sphinx-rtd-theme

**sphinx** packages may fail to install with pip, try with mamba::

    mamba install sphinx sphinx_rtd_theme


2. After commiting all changes, be sure that the release shows up on ``history.rst``
   and list all fixes and additions.

3. Update the version on ``locuaz/__init__.py``:

.. code-block:: console

    __version__ = "0.4.1"

and on ``setup.cfg``:

.. code-block:: console

    [metadata]
    name = locuaz
    version = 0.4.1
    author =  Patricio Barletta
    ...

4. Make sure you're on the right conda environment so you can build the documentation,
   the project and that ``twine`` doesn't report any errors on the distribution
   ``.tar.gaz`` nor the wheels::

    sphinx-build -b html docs/source/ docs/build/html/
    pyproject-build
    twine check dist/*

5. Commit ``locuaz/__init__.py``, and ``setup.cfg`` with a message along the lines of
   "Bump up version to X.X.X"::

    git add locuaz/__init__.py setup.cfg
    git commit -m "Bump up version to 0.4.1"
    git push origin main

6. Tag the release:

.. code-block:: console

    git tag 0.4.1
    git push origin --tags

7. The last step will trigger a GHAction to publish the project on pypi.
   Get the sha256 hash from the project's `pypi`_

8. Use the sha256 to update the version on the conda recipe:

.. code-block:: console

    {% set name = "locuaz" %}
    {% set version = "0.4.1" %}

    package:
      name: {{ name|lower }}
      version: {{ version }}

    source:
      url: https://pypi.io/packages/source/{{ name[0] }}/{{ name }}/locuaz-{{ version }}.tar.gz
      sha256: 93eed64eddba7ef6137d894aa26e3000ea51953153a26f2efeeb773990674f0a

9. Make sure you can build the conda recipe locally::

    conda mambabuild .

10. Make a PR to the feedstock and update it following conda-forge `instructions`_
11. Update the locuaz version on the apptainer definition file here::

        pip install locuaz==0.7.0 --root-user-action=ignore

And here::

    %labels
    Author Patricio Barletta
    Version 0.7.0

12. Build the container::

        sudo apptainer build locuaz.sif locuaz.def

13. Log into the GitHub container registry::

        apptainer remote login --username pgbarletta docker://ghcr.io

14. Upload the image to the GitHub registry::

        apptainer push locuaz.sif

Check `this`_ blog post for more info on apptainer.

.. _this: https://ana.run/blog/apptainer

Modifying the schema
---------------------
If at any time there's a change on the ``schema.yaml`` file, some documentation
needs to be updated:

1. ``configurationfile.rst``
2. Tutorials: ``tutorialsimple.rst``, ``tutorialtleap.rst`` and ``tutorialligand.rst``.
3. The configuration files in the example folders: ``/simple_tutorial/config_simple.yaml``,
   ``/tleap_tutorial/config_nb.yaml`` and ``/ligand_tutorial/config_ligand.yaml``.
   Then, update these new config files to the onedrive.

.. _pypi: https://pypi.org/project/locuaz/#files
.. _instructions: https://conda-forge.org/docs/maintainer/updating_pkgs.html

