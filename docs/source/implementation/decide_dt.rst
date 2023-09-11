
.. include:: /references.txt

.. _decide_dt:

##################
`src/decide_dt.c`_
##################

.. _src/decide_dt.c: https://github.com/NaokiHori/SimpleNSSolver/blob/main/src/decide_dt.c

This file contains functions which determine the next time step size ``dt`` (:math:`\Delta t`) used to integrate the governing equations.

.. mydeclare:: /../../src/decide_dt.c
   :language: c
   :tag: decide_dt

As discussed in :ref:`the temporal integration techniques <temporal_discretisation>`, there are two time scales in the current system which affect the stability to integrate the equations in time, which are considered separately here:

.. myliteralinclude:: /../../src/decide_dt.c
   :language: c
   :tag: compute advective and diffusive constraints

After all possible time step sizes (one scalar value from the advective terms, ``NDIMS`` values from the diffusive terms) are computed, I extract the smallest value

.. myliteralinclude:: /../../src/decide_dt.c
   :language: c
   :tag: choose smallest value as dt

and use it as the next time step.

.. note::

   When the diffusive terms are treated implicitly, I eliminate the constraint in the direction.
   Namely, I only take into account the diffusive constraints if the direction is treated explicitly.

*********************
Advective constraints
*********************

.. mydeclare:: /../../src/decide_dt.c
   :language: c
   :tag: decide_dt_adv

Advective restriction is computed here, in particular by computing the ratio of the local grid size and the magnitude of the local velocity:

* :math:`x` direction

   .. myliteralinclude:: /../../src/decide_dt.c
      :language: c
      :tag: compute grid-size over velocity in x

* :math:`y` direction

   .. myliteralinclude:: /../../src/decide_dt.c
      :language: c
      :tag: compute grid-size over velocity in y

* :math:`z` direction

   .. myliteralinclude:: /../../src/decide_dt.c
      :language: c
      :tag: compute grid-size over velocity in z

.. note::

   A small number is added to avoid zero divisions.

After unifying the result among all processes, a safety factor (given by the user) is multiplied:

.. myliteralinclude:: /../../src/decide_dt.c
   :language: c
   :tag: unify result, multiply safety factor

*********************
Diffusive constraints
*********************

.. mydeclare:: /../../src/decide_dt.c
   :language: c
   :tag: decide_dt_dif

Diffusive restrictions in all dimensions are computed here following

.. math::

   \frac{C}{2 \times \text{NDIMS}} \left[ \min \left( \Delta x_i \right) \right]^2,

in each dimension, where :math:`C` is a pre-factor (the inverse of the diffusivity times the safety factor) and :math:`\min \left( \Delta x_i \right)` is the smallest grid size in the direction:

.. myliteralinclude:: /../../src/decide_dt.c
   :language: c
   :tag: compute diffusive constraints

.. note::

   Roughly speaking the safety factors are smaller than 1.
   Although the three-step Runge-Kutta scheme which is used throughout this project allows slightly larger value than 1 (c.f. |COSTA2018|), one should remember that, the larger :math:`\Delta t` is used, the more dissipative the scheme becomes (c.f. |MORINISHI1998|).

