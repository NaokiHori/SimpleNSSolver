
.. _fluid:

############
`src/fluid`_
############

.. _src/fluid: https://github.com/NaokiHori/SimpleNSSolver/tree/main/src/fluid

This directory contains source files which are used to integrate the :ref:`governing equations <governing_equations>`.

.. toctree::
   :maxdepth: 1

   boundary
   compute_potential/main
   correct_velocity
   init
   integrate/main
   update_pressure

*******
fluid_t
*******

A structure ``fluid_t`` is defined in `include/fluid.h <https://github.com/NaokiHori/SimpleNSSolver/blob/main/include/fluid.h>`_, which contains all arrays and buffers which are used to integrate the momentum field:

.. myliteralinclude:: /../../include/fluid.h
   :language: c
   :tag: definition of a structure fluid_t

Each member plays the following role:

* ``u[xyz]``

   Velocity in the :math:`x`, :math:`y` and :math:`z` directions (:math:`\ux`, :math:`\uy` and :math:`\uz`)

* ``p``

   Pressure (:math:`p`, roughly speaking)

* ``psi``

   Scalar potential which projects the non-solenoidal velocity field to the solenoidal one (:math:`\psi`)

* ``t``

   Temperature field (or passive scalar field when the buoyancy force is neglected)

* ``src(u[xyz]|t)[abg]``

   Right-hand-side terms of the momentum and the internal-energy equations (i.e. source terms of the Runge-Kutta scheme).
   Suffices ``a``, ``b`` and ``g`` are used to distinguish source terms (implying :math:`\alpha`, :math:`\beta`, and :math:`\gamma`, respectively).

* ``m_dif``, ``t_dif``

   Since there are two fields, two parameters :math:`Ra` and :math:`Pr` are needed to identify the diffusive process.
   To avoid computing the diffusivities ``m_dif``:

   .. math::

      \frac{{\sqrt{Pr}}}{\sqrt{Ra}}

   and ``t_dif``:

   .. math::

      \frac{1}{\sqrt{Pr}} \frac{1}{\sqrt{Ra}}

   every time, they are pre-computed and stored:

   .. myliteralinclude:: /../../src/fluid/init.c
      :language: c
      :tag: compute diffusivities

.. seealso::

   :ref:`SMAC method <smac_method>`

   :ref:`Integrating the temperature field <temperature_integration>`

