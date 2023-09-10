
.. _fluid_update_pressure:

.. include:: /references.txt

##############################
`src/fluid/update_pressure.c`_
##############################

.. _src/fluid/update_pressure.c: https://github.com/NaokiHori/SimpleNSSolver/blob/main/src/fluid/update_pressure.c

Update the pressure field :math:`p` using the scalar potential :math:`\psi`.

.. mydeclare:: /../../src/fluid/update_pressure.c
   :language: c
   :tag: fluid_update_pressure

To update the pressure field, a scalar potential :math:`\psi`, which is computed by :ref:`solving the Poisson equation <fluid_compute_potential>`, is added:

.. math::

   p^{k+1}
   =
   p^k
   +
   \psi,

which is the *explicit* contribution.

.. myliteralinclude:: /../../src/fluid/update_pressure.c
   :language: c
   :tag: explicit contribution

As discussed in :ref:`the temporal discretisation <smac_method>`, when the diffusive terms are treated implicitly, additional terms appear, e.g.

.. math::

   p^{k+1}
   =
   p^k
   +
   \psi
   -
   \frac{\gamma^k \Delta t}{2} \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{\delta^2 \psi}{\delta x^2},

when the :math:`x` direction is treated implicitly (the *implicit* contributions).
The same applies to the other directions.

.. myliteralinclude:: /../../src/fluid/update_pressure.c
   :language: c
   :tag: x implicit contribution

.. myliteralinclude:: /../../src/fluid/update_pressure.c
   :language: c
   :tag: y implicit contribution

.. myliteralinclude:: /../../src/fluid/update_pressure.c
   :language: c
   :tag: z implicit contribution

The pre-factor

.. math::

   \frac{\gamma^k \Delta t}{2} \frac{\sqrt{Pr}}{\sqrt{Ra}}

is computed here:

.. myliteralinclude:: /../../src/fluid/update_pressure.c
   :language: c
   :tag: gamma dt diffusivity / 2

Also ``fluid->diffusivity`` is computed here:

.. myliteralinclude:: /../../src/fluid/init.c
   :language: c
   :tag: compute diffusivities

Finally the boundary condition is imposed on the new pressure field and the halo values are communicated:

.. myliteralinclude:: /../../src/fluid/update_pressure.c
   :language: c
   :tag: impose boundary conditions and communicate halo cells

.. note::

   In theory, the Helmholtz equations should be solved with respect to the new velocity field :math:`u_i^{k+1}`.
   For practical convenience, on the other hand, I solve them w.r.t. the predicted velocity field :math:`u_i^{*}`, which is the origin of the implicit contributions here.
   See :ref:`SMAC method <smac_method>` for more details.

   The implicit contributions, however, are known to play minor roles (e.g. |KAJISHIMA1999|).

