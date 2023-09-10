
.. _integrate:

##################
`src/integrate.c`_
##################

.. _src/integrate.c: https://github.com/NaokiHori/SimpleNSSolver/blob/main/src/integrate.c

``integrate.c`` contains a central function which integrates the :ref:`governing equations <governing_equations>`.

=======================
Deciding time step size
=======================

Before going into the Runge-Kutta iteration, the size of the time step ``dt`` is decided:

.. myliteralinclude:: /../../src/integrate.c
   :language: c
   :tag: decide time step size

======================
Runge-Kutta iterations
======================

The time integration is based on the three-step Runge-Kutta scheme.

.. note::

   At the beginning of each loop, I assume that the halo cells are communicated, and the proper boundary conditions are already enforced.

   To achieve the second-order accuracy in time, the following process is repeated for three times to update the scalar fields from step :math:`n` to step :math:`n+1`.

#. Update velocity field

   The velocities and the temperature field are updated by integrating the momentum and the internal-energy equations:

   .. myliteralinclude:: /../../src/integrate.c
      :language: c
      :tag: predict flow field

   Since I weakly couple the temperature and the momentum fields, the previous temperature field is used to compute the buoyancy force, which is added as a body force in the momentum equation.

   In the first part :ref:`fluid_compute_rhs <fluid_compute_rhs>`, the right-hand-side terms of the equations are computed, which are used to actually update the velocity field in the second part :ref:`fluid_predict_field <fluid_predict_field>`.

   .. seealso::

      :ref:`Numerical method <numerics>` for more details.

#. Compute scalar potential and correct velocity

   Usually, the new velocity field computed in the previous step does not satisfy the mass conservation.
   To enforce the incompressibility, a scalar potential :math:`\psi` is computed by solving a Poisson equation here:

   .. myliteralinclude:: /../../src/integrate.c
      :language: c
      :tag: compute scalar potential

   which is used to correct the velocity field to be solenoidal:

   .. myliteralinclude:: /../../src/integrate.c
      :language: c
      :tag: correct velocity field to satisfy mass conservation

#. Update pressure

   Finally the pressure is updated:

   .. myliteralinclude:: /../../src/integrate.c
      :language: c
      :tag: update pressure

.. seealso::

   :ref:`SMAC method <smac_method>`

   :ref:`Integrating the temperature field <temperature_integration>`

