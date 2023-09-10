
.. _fluid_boundary:

#####################
`src/fluid/boundary`_
#####################

.. _src/fluid/boundary: https://github.com/NaokiHori/SimpleNSSolver/tree/main/src/fluid/boundary

This directory contains functions to update the boundary values of the velocity, the pressure, the scalar potential, and the temperature.
In particular, there are two major roles which are processed here:

* imposing given boundary conditions,

* exchanging the halo cells using the inter-process communications.

.. mydeclare:: /../../src/fluid/boundary/ux.c
   :language: c
   :tag: fluid_update_boundaries_ux

.. mydeclare:: /../../src/fluid/boundary/uy.c
   :language: c
   :tag: fluid_update_boundaries_uy

.. mydeclare:: /../../src/fluid/boundary/uz.c
   :language: c
   :tag: fluid_update_boundaries_uz

.. mydeclare:: /../../src/fluid/boundary/p.c
   :language: c
   :tag: fluid_update_boundaries_p

.. mydeclare:: /../../src/fluid/boundary/psi.c
   :language: c
   :tag: fluid_update_boundaries_psi

.. mydeclare:: /../../src/fluid/boundary/t.c
   :language: c
   :tag: fluid_update_boundaries_t

*******************
Boundary conditions
*******************

Since the domain is wall-bounded in :math:`x`,

* ``ux``: impermeable condition,

* ``uy``: no-slip condition,

* ``uz``: no-slip condition,

* ``p``: Neumann condition,

* ``psi``: Neumann condition,

* ``t``: Dirichlet condition,

are to be imposed on the walls:

.. myliteralinclude:: /../../src/fluid/boundary/ux.c
   :language: c
   :tag: set boundary values

.. myliteralinclude:: /../../src/fluid/boundary/uy.c
   :language: c
   :tag: set boundary values

.. myliteralinclude:: /../../src/fluid/boundary/uz.c
   :language: c
   :tag: set boundary values

.. myliteralinclude:: /../../src/fluid/boundary/p.c
   :language: c
   :tag: set boundary values

.. myliteralinclude:: /../../src/fluid/boundary/psi.c
   :language: c
   :tag: set boundary values

.. myliteralinclude:: /../../src/fluid/boundary/t.c
   :language: c
   :tag: set boundary values

*******************
Halo communications
*******************

Finite-difference methods evaluate the derivatives using the surrounding values.
Since I decompose the domain into several chunks, to evaluate the derivatives close to the edges, I need information which the process does not have by default.

Basically, when updating the momentum or the temperature fields, each process is responsible for updating the scalar values between ``1`` and ``jsize`` in the :math:`y` direction (see :ref:`src/fluid/predict <fluid_compute_rhs>`).
Thus communications are needed to obtain values at ``j = 0`` and ``jsize+1``, which are necessary to evaluate derivatives at ``j = 1`` and ``jsize``, respectively.

.. seealso::

   :ref:`src/halo.c <halo>`

