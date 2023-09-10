
.. _governing_equations:

#########
Equations
#########

*******************
Governing equations
*******************

================
Dimensional form
================

In this project, I consider the conservation laws of the mass, the momentum, and the internal energy, which are governed by

.. math::

   \der{\tilde{u}_i}{\tilde{x}_i}
   =
   0,

.. math::

   \der{\tilde{u}_i}{\tilde{t}}
   +
   \tilde{u}_j \der{\tilde{u}_i}{\tilde{x}_j}
   =
   -
   \frac{1}{\tilde{\rho}} \der{\tilde{p}}{\tilde{x}_i}
   +
   \tilde{\nu} \der{}{\tilde{x}_j} \der{\tilde{u}_i}{\tilde{x}_j}
   +
   \tilde{g}_i,

and

.. math::

   \der{\tilde{T}}{\tilde{t}}
   +
   \tilde{u}_i \der{\tilde{T}}{\tilde{x}_i}
   =
   \tilde{\kappa} \der{}{\tilde{x}_i} \der{\tilde{T}}{\tilde{x}_i},

respectively.

.. note::

   Physical properties (e.g. the density :math:`\tilde{\rho}`, the kinematic viscosity :math:`\tilde{\nu}`, the thermal diffusivity :math:`\tilde{\kappa}`) are assumed to be constant in time and in space.

   :math:`\tilde{q}` implies that the quantity :math:`q` is dimensional (i.e. before normalised).

Also, by combining the mass balance and the momentum balance, the equation of the squared velocity

.. math::

   \der{\tilde{k}}{\tilde{t}}
   +
   \tilde{u}_i \der{\tilde{k}}{\tilde{x}_i}
   =
   -
   \tilde{u}_i \frac{1}{\tilde{\rho}} \der{\tilde{p}}{\tilde{x}_i}
   +
   \tilde{\nu} \der{}{\tilde{x}_j} \left( \tilde{u}_i \der{\tilde{u}_i}{\tilde{x}_j} \right)
   -
   \tilde{\nu} \der{\tilde{u}_i}{\tilde{x}_j} \der{\tilde{u}_i}{\tilde{x}_j}
   +
   \tilde{u}_i \tilde{g}_i

comes out, where

.. math::

   \tilde{k}
   \equiv
   \frac{1}{2}
   \tilde{u}_i \tilde{u}_i.

Similarly I obtain the equation of the squared temperature

.. math::

   \der{\tilde{h}}{\tilde{t}}
   +
   \tilde{u}_i \der{\tilde{h}}{\tilde{x}_i}
   =
   \tilde{\kappa} \der{}{\tilde{x}_i} \left( \tilde{T} \der{\tilde{T}}{\tilde{x}_i} \right)
   -
   \tilde{\kappa} \der{\tilde{T}}{\tilde{x}_i} \der{\tilde{T}}{\tilde{x}_i}

using the internal energy balance, where

.. math::

   \tilde{h}
   \equiv
   \frac{1}{2}
   \tilde{T} \tilde{T}.

These additional equations also play important roles in the following discussion.

Although the derivations can be found everywhere, they are briefly discussed here for the sake of completeness.

.. toctree::
   :maxdepth: 1

   eq/mass
   eq/momentum
   eq/internal_energy

====================
Non-dimensional form
====================

In this project, I focus on `Rayleigh-Bénard convection <https://en.wikipedia.org/wiki/Rayleigh–Bénard_convection>`_, which is an excellent model problem to shed light on the conservation properties.
By adopting `Boussinesq approximation <https://en.wikipedia.org/wiki/Boussinesq_approximation_(buoyancy)>`_ and normalise the equations with proper scales, I obtain the following non-dimensional equations which play the central role in this project.

.. _eq_mass:

* Mass balance

   .. math::

      \der{u_i}{x_i}
      =
      0.

.. _eq_momentum:

* Momentum balance

   .. math::

      \der{u_i}{t}
      +
      u_j \der{u_i}{x_j}
      =
      -
      \der{p}{x_i}
      +
      \frac{\sqrt{Pr}}{\sqrt{Ra}} \der{}{x_j} \der{u_i}{x_j}
      +
      T \delta_{ix}.

.. _eq_temperature:

* Internal energy balance

   .. math::

      \der{T}{t}
      +
      u_i \der{T}{x_i}
      =
      \frac{1}{\sqrt{Pr} \sqrt{Ra}} \der{}{x_i} \der{T}{x_i}.

.. _eq_squared_velocity:

* Equation of the squared velocity

   .. math::

      \der{k}{t} + u_i \der{k}{x_i}
      =
      -
      u_i \der{p}{x_i}
      +
      \frac{\sqrt{Pr}}{\sqrt{Ra}} \der{}{x_j} \left( u_i \der{u_i}{x_j} \right)
      -
      \frac{\sqrt{Pr}}{\sqrt{Ra}} \der{u_i}{x_j} \der{u_i}{x_j}
      +
      u_i T \delta_{ix}.

.. _eq_squared_temperature:

* Equation of the squared temperature

   .. math::

      \der{h}{t}
      +
      u_i \der{h}{x_i}
      =
      \frac{1}{\sqrt{Pr} \sqrt{Ra}} \der{}{x_i} \left( T \der{T}{x_i} \right)
      -
      \frac{1}{\sqrt{Pr} \sqrt{Ra}} \der{T}{x_i} \der{T}{x_i}.

Here Rayleigh number :math:`Ra` and Prandtl number :math:`Pr` are dimensionless parameters given by

.. math::

   Ra & = \frac{\tilde{\beta} \tilde{g} {\tilde{l_x}}^3 \left( \Delta \tilde{T} \right)}{\tilde{\nu} \tilde{\kappa}}, \\
   Pr & = \frac{\tilde{\nu}}{\tilde{\kappa}},

where :math:`\tilde{\beta}`, :math:`\tilde{g}`, :math:`\tilde{l_x}`, and :math:`\Delta \tilde{T} = \tilde{T}_{H} - \tilde{T}_{L}` are the thermal expansion coefficient :math:`\left[ K^{-1} \right]`, the gravitational acceleration :math:`\left[ L T^{-2} \right]`, the distance between the walls :math:`\left[ L \right]`, and the temperature difference :math:`\left[ K \right]`, respectively.

.. image:: image/schematic.png
   :align: center
   :width: 400

Periodic boundary conditions are imposed in the homogeneous directions :math:`y` and :math:`z`.
The boundary conditions in the :math:`x` (wall-normal) direction are listed here:

* :math:`\ux = 0`: Dirichlet condition, impermeable walls.

* :math:`\uy = \uz = 0`: Dirichlet condition, no-slip and stationary walls.

* :math:`\partial p / \partial x = 0`: Neumann condition.

* :math:`\vat{T}{x = 0}, \vat{T}{x = l_x \equiv 1}`: Dirichlet condition, fixed temperature satisfying :math:`\vat{T}{x = 0} - \vat{T}{x = 1} = 1`.

.. note::

   * Without loss of generality, I can fix :math:`\tilde{\beta}`, :math:`\tilde{g}`, :math:`\tilde{l_x}`, and :math:`\Delta \tilde{T}` to unity, which is assumed in this project.

   * The reference velocity scale :math:`\tilde{U} \left[ L T^{-1} \right]` is defined by the other parameters :math:`\tilde{U} = \sqrt{\tilde{\beta} \tilde{g} \tilde{l_x} \left( \Delta \tilde{T} \right)} \left( = 1 \right)`, which is often called as the free-fall velocity.

   * Although it is numerically trivial to let the walls move in the homogeneous directions, the following discussion assumes the walls are at rest.
     In particular, I need to modify :ref:`the Nusselt number based on the kinetic energy dissipation <nu_kinetic_energy_dissipation>` if the walls move.

**************
Nusselt number
**************

One of the most important features of the Rayleigh-Bénard convections is the theoretical relationships of the Nusselt number :math:`Nu`, which measures how much the heat transfer between the two walls is enhanced by the convective effects.
In this part, I mathematically define the heat flux and the Nusselt number, which is followed by the derivations of several :math:`Nu` relationships in the continuous domain.

.. toctree::
   :maxdepth: 1

   nu/heat_flux
   nu/kinetic_energy_injection
   nu/kinetic_energy_dissipation
   nu/thermal_energy_dissipation

.. seealso::

   Although all :math:`Nu` should give the identical result in the continuous domain, it is non-trivial to satisfy these relations umerically (discrete domain).
   This is discussed in detail in :ref:`the numerical method <numerics>`.

