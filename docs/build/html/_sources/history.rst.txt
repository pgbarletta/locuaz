=======
History
=======

0.3.8 (2023-05-)
------------------
 * *DLPacker* data files ``library.npz`` and ``charges.rtp`` are now included with the install. Only the weights have
   to be downloaded and extracted into a dir whose path must be specified in the ``config['paths']['mutator']`` option.

0.3.7 (2023-05-12)
------------------
 * Added Directed Acyclic Graph tracking of the protocol, so a plot of the progression of the protocol can be done,
   both of the iteration names and the mutations performed on each mutation.
 * Added docs on https://locuaz.readthedocs.io/
 * Made DLPacker part of the repo. Used for performing mutations.
 * Added metropolis Pruner.

0.2.1 (2023-04-20)
------------------
* The protocol is now fully installable by pip, provided that ambertools and tensorflow are present in the conda environment (no available pip install for them)


0.2.0 (2022-05-13)
------------------
* First fully functional release.

0.1.0 (2022-05-25)
------------------
* First release on PyPI.
