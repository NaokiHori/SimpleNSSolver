########
Momentum
########

.. mydeclare:: /../../src/logging/momentum.c
   :language: c
   :tag: logging_check_momentum

This function computes the net momentum in each direction (e.g. :math:`\int \ux dx dy dz`) and writes them to a file.

The discrete form leads to

.. math::

   \sum_{i} \sum_{j} \sum_{k}
   \vat{
      \left(
         \ux
         \Delta x
         \Delta y
         \Delta z
      \right)
   }{\xic, \xjc, \xkc},

which are implemented as

.. myliteralinclude:: /../../src/logging/momentum.c
   :language: c
   :tag: compute total x-momentum

.. myliteralinclude:: /../../src/logging/momentum.c
   :language: c
   :tag: compute total y-momentum

.. myliteralinclude:: /../../src/logging/momentum.c
   :language: c
   :tag: compute total z-momentum

