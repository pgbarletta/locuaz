Pruners
====================
This block is in charge of deciding which *Iterations* go through the next epoch, according to their scores.

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
For a new complex to pass onto the next *epoch*, it has to beat all previous complexes that gave rise to it.
Using this pruner, a new complex beats a previous one when the means of its scores are lower than those
from the previous complex, given that a lower score number indicates higher affinity between the
target and the binder.

Now, since scoring functions tend to diverge, it's not necessary for all of them to improve.
The user can set a minimum threshold of scoring functions that have to improve so a binder can
be considered better than another.

For example, if a user is using **4** scoring functions, the user can set the ``threshold``
to, say, **2** and if 2 scoring functions indicate an improvement, then the new complex
beats the old one. More formally:

.. math::

   c_{k}^{i} = \left\{
    \begin{array}{l}
    1, \quad \text{if} \ avg\left(score^{i+1}_{k}\right) < avg\left(score^{i}_{k}\right) \\
    0, \quad \text{else}
    \end{array}
    \right.

Where
:math:`c_{k}^{i}` counts how many scoring functions improved with respect to the previous complex,
while :math:`k` and :math:`i` are the indices for scoring functions and complexes, respectively.


Then, a consensus number :math:`C^{i}` is obtained by adding all the :math:`c_{k}^{i}` for the
:math:`N` scoring functions, and this number is compared against the threshold :math:`T`:

.. math::

    C^{i} = \sum_{k=1}^{N} c_{k}^{i}
    
.. math::
    C^{i} \geq T ?

If the last statement is true, then complex :math:`i+1` beats :math:`i` and can be considered
for a next round of mutations.

.. automodule:: locuaz.prunerconsensus
   :members:
   :undoc-members:
   :show-inheritance:

locuaz.prunermetropolis module
--------------------------------
When using only 1 scoring functino, the well known *metropolis acceptance criteria* can be used to decide
whether a complex passes to the next *epoch*:

.. math::

    \text{Acceptance ratio} = \min\left(1, \exp\left(-\frac{\Delta E}{k_B T}\right)\right)

Then a random number between :math:`0` and :math:`1` is generated and if the *Acceptance ratio* is above it,
then the new complex is considered to beat the old one.

.. automodule:: locuaz.prunermetropolis
   :members:
   :undoc-members:
   :show-inheritance:

