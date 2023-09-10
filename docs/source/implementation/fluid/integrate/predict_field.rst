
.. _fluid_predict_field:

######################################
`src/fluid/integrate/predict_field.c`_
######################################

.. _src/fluid/integrate/compute_rhs.c: https://github.com/NaokiHori/SimpleNSSolver/blob/main/src/fluid/integrate/predict_field.c

Functions to update the momentum field

.. math::

   u_i^n \rightarrow u_i^*

and the temperature field

.. math::

   T^n \rightarrow T^{n+1}

are discussed here.

See :ref:`the temporal discretisation <temporal_discretisation>` and :ref:`the spatial discretisation <spatial_discretisation>` for details.

.. note::

   See :ref:`fluid/update_pressure.c <fluid_update_pressure>` for the pressure field, which is not discussed here since the pressure is assumed to be a scalar potential to enforce the incompressibility.

.. mydeclare:: /../../src/fluid/integrate/predict_field.c
   :language: c
   :tag: fluid_predict_field

The right-hand-side terms of the Runge-Kutta scheme are already computed by :ref:`fluid_compute_rhs <fluid_compute_rhs>`.
Here, those terms are adopted to integrate the scalar fields in time.

.. note::

   Although :math:`T` is used throughout this page for notational simplicity, the same applies to the other quantities (i.e. velocity in each direction).

First of all, the increment of the scalar field :math:`\delta T \equiv T^{k+1} - T^{k}` is computed:

.. math::

   \delta T = \left(
        \alpha^k \mathcal{A}
      + \beta^k  \mathcal{B}
      + \gamma^k \mathcal{G}
   \right) \Delta t,

where :math:`\mathcal{A}`, :math:`\mathcal{B}`, and :math:`\mathcal{G}` are source terms which are already computed in ``compute_src`` and correspond to ``srca``, ``srcb``, and ``srcg``, respectively.

As discussed in :ref:`the temporal discretisation <temperature_integration>`, the update process differs for different diffusive treatments.
When all diffusive terms are treated explicitly in time, I simply have

.. math::

   T^{k+1} = T^k + \delta T.

When the diffusive terms are partially treated implicitly (e.g. only :math:`x` direction), I have

.. math::

   T^{k+1} = T^k + \left( 1 - \frac{\gamma^k \Delta t}{2 \sqrt{Pr} \sqrt{Ra}} \frac{\delta^2}{\delta x^2} \right)^{-1} \delta T,

where linear systems in :math:`x` direction needs to be solved before being added to :math:`T^k`.

When all diffusive terms are treated implicitly, I have

.. math::

   T^{k+1}
   =
   T^k
   +
   \left( 1 - \frac{\gamma^k \Delta t}{2 \sqrt{Pr} \sqrt{Ra}} \frac{\delta^2}{\delta y^2} \right)^{-1}
   \left( 1 - \frac{\gamma^k \Delta t}{2 \sqrt{Pr} \sqrt{Ra}} \frac{\delta^2}{\delta x^2} \right)^{-1}
   \delta T,

or

.. math::

   T^{k+1}
   =
   T^k
   +
   \left( 1 - \frac{\gamma^k \Delta t}{2 \sqrt{Pr} \sqrt{Ra}} \frac{\delta^2}{\delta z^2} \right)^{-1}
   \left( 1 - \frac{\gamma^k \Delta t}{2 \sqrt{Pr} \sqrt{Ra}} \frac{\delta^2}{\delta y^2} \right)^{-1}
   \left( 1 - \frac{\gamma^k \Delta t}{2 \sqrt{Pr} \sqrt{Ra}} \frac{\delta^2}{\delta x^2} \right)^{-1}
   \delta T,

where linear systems in each direction should be solved before being added to :math:`T^k`.

On the basis of these relationships, the procedure for updating the field is as follows.

#. Compute :math:`\delta T`

   The temperature increment :math:`\delta T` is computed and assigned to a buffer:

   .. myliteralinclude:: /../../src/fluid/integrate/ux.c
      :language: c
      :tag: compute increments

   .. myliteralinclude:: /../../src/fluid/integrate/uy.c
      :language: c
      :tag: compute increments

   .. myliteralinclude:: /../../src/fluid/integrate/uz.c
      :language: c
      :tag: compute increments

   .. myliteralinclude:: /../../src/fluid/integrate/t.c
      :language: c
      :tag: compute increments

#. Solve linear systems in each direction (**only when** the diffusive term in the direction is treated implicitly)

   See :ref:`src/linear_system.c <linear_system>`.

#. Update :math:`T^k` to :math:`T^{k+1}`

   Now the increment :math:`\delta T` is ready, which is added to the temperature ``temp``:

   .. myliteralinclude:: /../../src/fluid/integrate/ux.c
      :language: c
      :tag: the field is actually updated

   .. myliteralinclude:: /../../src/fluid/integrate/uy.c
      :language: c
      :tag: the field is actually updated

   .. myliteralinclude:: /../../src/fluid/integrate/uz.c
      :language: c
      :tag: the field is actually updated

   .. myliteralinclude:: /../../src/fluid/integrate/t.c
      :language: c
      :tag: the field is actually updated

