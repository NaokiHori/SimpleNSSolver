
.. _fluid_correct_velocity:

#############################
`src/fluid/correct_velocity`_
#############################

.. _src/fluid/correct_velocity: https://github.com/NaokiHori/SimpleNSSolver/tree/main/src/fluid/correct_velocity

This directory contains functions to project (correct) the non-solenoidal velocity :math:`u_i^*` to the solenoidal (and new-step) field :math:`u_i^{n+1}` following the relation

.. math::

   u_i^{n+1}
   -
   u_i^*
   =
   - \frac{\delta \psi}{\delta x_i} \gamma \Delta t,

where :math:`\psi` is a scalar potential which was computed by :ref:`fluid_compute_potential <fluid_compute_potential>`.

``fluid_correct_velocity`` is a main function, which calls sub functions which update each velocity field.

.. seealso::

   :ref:`SMAC method <smac_method>`.

.. mydeclare:: /../../src/fluid/correct_velocity/ux.c
   :language: c
   :tag: fluid_correct_velocity_ux

.. mydeclare:: /../../src/fluid/correct_velocity/uy.c
   :language: c
   :tag: fluid_correct_velocity_uy

.. mydeclare:: /../../src/fluid/correct_velocity/uz.c
   :language: c
   :tag: fluid_correct_velocity_uz

The pre-factor is computed here:

.. myliteralinclude:: /../../src/fluid/correct_velocity/main.c
   :language: c
   :tag: compute prefactor gamma dt

The above equation is directly discretised and implemented:

.. myliteralinclude:: /../../src/fluid/correct_velocity/ux.c
   :language: c
   :tag: correct x velocity

.. myliteralinclude:: /../../src/fluid/correct_velocity/uy.c
   :language: c
   :tag: correct y velocity

.. myliteralinclude:: /../../src/fluid/correct_velocity/uz.c
   :language: c
   :tag: correct z velocity

Also the boundary values are imposed and the halo cells are updated:

.. myliteralinclude:: /../../src/fluid/correct_velocity/ux.c
   :language: c
   :tag: update boundary and halo cells

.. myliteralinclude:: /../../src/fluid/correct_velocity/uy.c
   :language: c
   :tag: update boundary and halo cells

.. myliteralinclude:: /../../src/fluid/correct_velocity/uz.c
   :language: c
   :tag: update boundary and halo cells

.. seealso::

   :ref:`Boudary treatments <fluid_boundary>`.

