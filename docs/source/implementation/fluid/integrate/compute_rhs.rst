
.. _fluid_compute_rhs:

####################################
`src/fluid/integrate/compute_rhs.c`_
####################################

.. _src/fluid/integrate/compute_rhs.c: https://github.com/NaokiHori/SimpleNSSolver/blob/main/src/fluid/integrate/compute_rhs.c

Functions to compute the right-hand-side terms of the momentum and the internal energy equation are described here.

See :ref:`the temporal discretisation <temporal_discretisation>` and :ref:`the spatial discretisation <spatial_discretisation>` for details.

.. mydeclare:: /../../src/fluid/integrate/compute_rhs.c
   :language: c
   :tag: fluid_compute_rhs

To update a scalar field from :math:`\left\{ \cdots \right\}^k` to :math:`\left\{ \cdots \right\}^{k+1}` following the method discussed in :ref:`the temporal discretisation <temperature_integration>`, I need the right-hand-side terms: the advective terms :math:`A^k`, and the diffusive terms :math:`D^k`.
I also need the pressure-gradient contribution :math:`P^k` to the momentum balance, and in addition the buoyancy forcing is coupled in the :math:`x` direction.

The following three buffers are used:

.. list-table:: Buffers
   :widths: 25 25 50
   :header-rows: 1

   * - Name
     - Pre-factor
     - Description
   * - ``srca``
     - :math:`\alpha^k`
     - Explicit terms, current  contribution
   * - ``srcb``
     - :math:`\beta^k`
     - Explicit terms, previous contribution
   * - ``srcg``
     - :math:`\gamma^k`
     - Implicit terms

First of all,

   * ``srca``, which contains the information of the previous Runge-Kutta step, is copied to ``srcb``, so that the new values can be assigned to ``srca``.

   * ``srca`` and ``srcg`` are zero-cleared.

.. myliteralinclude:: /../../src/fluid/integrate/compute_rhs.c
   :language: c
   :tag: copy previous k-step source term and reset

Note that the zero-clearing is not needed for :math:`\beta` buffers even at :math:`k = 0` since the Runge-Kutta coefficient :math:`\beta^0` is 0.

Then all different advective and diffusive contributions are added to those buffers.

.. note::

   The mathematical forms and their derivations are extensively discussed in :ref:`the spatial discretisation <spatial_discretisation>`.
   In short, the schemes are designed to keep the properties which the original governing equations have.

   Hereafter, superscripts :math:`k`, which denote the Runge-Kutta iteration and should be on all scalar fields, are dropped for notational simplicity.

* Advective contributions

   Since they are always treated explicitly in time, they are added to ``srca``.

   .. myliteralinclude:: /../../src/fluid/integrate/ux.c
      :language: c
      :tag: advective contributions, always explicit

   * :ref:`Advection of x momentum in x <impl_adv_x_x>`
   * :ref:`Advection of x momentum in y <impl_adv_x_y>`
   * :ref:`Advection of x momentum in z <impl_adv_x_z>`

   .. myliteralinclude:: /../../src/fluid/integrate/uy.c
      :language: c
      :tag: advective contributions, always explicit

   * :ref:`Advection of y momentum in x <impl_adv_y_x>`
   * :ref:`Advection of y momentum in y <impl_adv_y_y>`
   * :ref:`Advection of y momentum in z <impl_adv_y_z>`

   .. myliteralinclude:: /../../src/fluid/integrate/uz.c
      :language: c
      :tag: advective contributions, always explicit

   * :ref:`Advection of z momentum in x <impl_adv_z_x>`
   * :ref:`Advection of z momentum in y <impl_adv_z_y>`
   * :ref:`Advection of z momentum in z <impl_adv_z_z>`

   .. myliteralinclude:: /../../src/fluid/integrate/t.c
      :language: c
      :tag: advective contributions, always explicit

   * :ref:`Advection of T in x <impl_adv_t_x>`
   * :ref:`Advection of T in y <impl_adv_t_y>`
   * :ref:`Advection of T in z <impl_adv_t_z>`

* Diffusive contributions

   The treatment is user-dependent and determined at runtime (see ``exec.sh``).
   When the term is treated explicitly in time, it is added to ``srca``; otherwise to ``srcg``:

   .. myliteralinclude:: /../../src/fluid/integrate/ux.c
      :language: c
      :tag: diffusive contributions, can be explicit or implicit

   * :ref:`Diffusion of x momentum in x <impl_dif_x_x>`
   * :ref:`Diffusion of x momentum in y <impl_dif_x_y>`
   * :ref:`Diffusion of x momentum in z <impl_dif_x_z>`

   .. myliteralinclude:: /../../src/fluid/integrate/uy.c
      :language: c
      :tag: diffusive contributions, can be explicit or implicit

   * :ref:`Diffusion of y momentum in x <impl_dif_y_x>`
   * :ref:`Diffusion of y momentum in y <impl_dif_y_y>`
   * :ref:`Diffusion of y momentum in z <impl_dif_y_z>`

   .. myliteralinclude:: /../../src/fluid/integrate/uz.c
      :language: c
      :tag: diffusive contributions, can be explicit or implicit

   * :ref:`Diffusion of z momentum in x <impl_dif_z_x>`
   * :ref:`Diffusion of z momentum in y <impl_dif_z_y>`
   * :ref:`Diffusion of z momentum in z <impl_dif_z_z>`

   .. myliteralinclude:: /../../src/fluid/integrate/t.c
      :language: c
      :tag: diffusive contributions, can be explicit or implicit

   * :ref:`Diffusion of T in x <impl_dif_t_x>`
   * :ref:`Diffusion of T in y <impl_dif_t_y>`
   * :ref:`Diffusion of T in z <impl_dif_t_z>`

   All diffusive terms contain the discrete Laplace operators, which are constant in time.
   Thus I compute them when the simulation is launched (see :ref:`domain/ <domain>`) and reuse the results.

   Also the diffusivities for the momentum equation

   .. math::

      \frac{\sqrt{Pr}}{\sqrt{Ra}}

   and for the internal energy

   .. math::

      \frac{1}{\sqrt{Pr}\sqrt{Ra}}

   are constant in time (see :ref:`the governing equations <governing_equations>`), which is pre-computed:

   .. myliteralinclude:: /../../src/fluid/init.c
      :language: c
      :tag: compute diffusivities

* Buoyancy (:math:`x` momentum equation)

   The behaviour of this part is dominated by a flag ``param_add_buoyancy`` (see :ref:`src/param <param>`), which determines whether the buoyancy force under the Boussinesq approximation is added to the right-hand side of the :math:`x` momentum equation (``true``) or not (``false``).
   When ``param_add_buoyancy`` is ``false``, the temperature field behaves as a passive scalar field.

   When ``param_add_buoyancy`` is ``true``, as discussed in :ref:`the governing equations <governing_equations>`, the buoyancy force results in :math:`T` (a simple body force term) in the :math:`x` momentum equation.
   Since :math:`T` are defined at cell centers, I need to interpolate them to the :math:`x` cell faces (where :math:`\ux` is located), where arithmetic average should be used.

   .. note::

      One might be tempted to use values which are linearly interpolated or weighted by cell volume, which breaks the energy balance.
      See :ref:`the spatial discretisation <spatial_discretisation>` for details, in particular :ref:`Nusselt number relations <nusselt_number_relations>`.

