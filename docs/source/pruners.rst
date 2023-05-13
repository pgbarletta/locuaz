Pruners
====================
These are the units in charge of deciding which Iterations go through the next epoch, according to their scores.

locuaz.pruners module
---------------------

.. automodule:: locuaz.pruners
   :members:
   :undoc-members:
   :show-inheritance:


locuaz.abstractpruner module
-----------------------------

.. automodule:: locuaz.abstractpruner
   :members:
   :undoc-members:
   :show-inheritance:

locuaz.prunerconsensus module
------------------------------

as

.. math::

   c_{k}^{i} = \left\{
    \begin{array}{l}
    1, \quad \text{if} \ avg\left(score^{i+1}_{k}\right) > avg\left(score^{i}_{k}\right) \\
    0, \quad \text{else}
    \end{array}
    \right.

Decide

.. math::

    C^{i} = \sum_{k=1}^{N} c_{k}^{i} \geq T

.. automodule:: locuaz.prunerconsensus
   :members:
   :undoc-members:
   :show-inheritance:

locuaz.prunermetropolis module
--------------------------------

.. automodule:: locuaz.prunermetropolis
   :members:
   :undoc-members:
   :show-inheritance:

