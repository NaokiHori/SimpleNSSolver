############################
Derivation: momentum balance
############################

.. include:: tilde_note.rst

********
Momentum
********

The momentum of the liquid in a control volume leads

.. math::

   \int_{V} \rho u_i dV \,\, \left[ M L T^{-1} \right].

The balance of this quantity is written as

.. math::

   \der{}{t} \int_{V} \rho u_i dV
   =
   -
   \int_{V} \der{u_j \rho u_i}{x_j} dV
   +
   \int_{\partial V} \sigma_{ij} n_j dS
   +
   \int_{V} \rho g_i dV, \,\, \left[ M L T^{-2} \right]

indicating that changes in the momentum are caused by

* the convection,

* the normal and the tangential surface forces,

* the body force.

Using the divergence theorem, the surface force leads to

.. math::

   \int_{V} \der{\sigma_{ij}}{x_j} dV,

and thus

.. math::

   \der{\rho u_i}{t}
   +
   \der{u_j \rho u_i}{x_j}
   =
   \der{\sigma_{ij}}{x_j}
   +
   \rho g_i,

where I assume the temporal derivative and the spatial integrals are commutative, and request the integrand to be zero.

Here, :math:`\sigma_{ij}` is the Cauchy stress tensor, which is regarded as

.. math::

   \sigma_{ij}
   =
   -
   p \delta_{ij}
   +
   2 \mu e_{ij},

where the incompressibility

.. math::

   \der{u_i}{x_i}
   =
   0

is assumed and :math:`\mu` is the dynamic viscosity of the liquid.
Also, :math:`p` is the reduced pressure and :math:`e_{ij}` is the symmetric part of the strain-rate tensor:

.. math::

   e_{ij}
   =
   \frac{1}{2} \left( \der{u_i}{x_j} + \der{u_j}{x_i} \right)

to close the equation.

Note that, using the incompressibility, I can write the advective term as

.. math::

   u_j \der{\rho u_i}{x_j}.

Finally, assuming that the physical properties are constant, I obtain

.. math::

   \der{u_i}{t}
   +
   u_j \der{u_i}{x_j}
   =
   -
   \frac{1}{\rho} \der{p}{x_i}
   +
   \nu \der{}{x_j} \der{u_i}{x_j}
   +
   g_i,

where :math:`\nu` is the kinematic viscosity defined as :math:`\mu / \rho`.

.. seealso::

   * `Navier-Stokes equations <https://en.wikipedia.org/wiki/Navier–Stokes_equations>`_
   * `Newtonian fluid <https://en.wikipedia.org/wiki/Newtonian_fluid>`_

****************
Squared velocity
****************

I consider to take the inner product of the velocity vector :math:`u_i` and the momentum balance derived above.

* Temporal derivative

   .. math::

      u_i \der{u_i}{t}
      =
      \der{u_i u_i}{t}
      -
      \der{u_i}{t} u_i,

   giving

   .. math::

      u_i \der{u_i}{t}
      =
      \frac{1}{2} \der{u_i u_i}{t}.

* Advective terms

   .. math::

      u_i u_j \der{u_i}{x_j}
      =
      u_j \der{u_i u_i}{x_j}
      -
      u_j \der{u_i}{t} u_i,

   giving

   .. math::

      u_i u_j \der{u_i}{x_j}
      =
      \frac{1}{2} u_j \der{u_i u_i}{x_j}.

* Pressure-gradient term

   .. math::

      -
      u_i \frac{1}{\rho} \der{p}{x_i}.

* Diffusive terms

   .. math::

      u_i \nu \der{}{x_j} \der{u_i}{x_j}
      =
      \nu \der{}{x_j} \left( u_i \der{u_i}{x_j} \right)
      -
      \nu \der{u_i}{x_j} \der{u_i}{x_j}.

   .. note::

      The first term in the right-hand side yields

      .. math::

         \nu \der{}{x_j} \left\{ \der{}{x_j} \left( u_i u_i \right) \right\}
         -
         \nu \der{}{x_j} \left( \der{u_i}{x_j} u_i \right),

      and thus I notice

      .. math::

         \nu \der{}{x_j} \left\{ \der{}{x_j} \left( \frac{1}{2} u_i u_i \right) \right\},

      which describes the diffusive effect of the squared velocity.

* Body force term

   .. math::

      u_i g_i.

Thus, I have

.. math::

   \der{k}{t}
   +
   u_i \der{k}{x_i}
   =
   -
   u_i \frac{1}{\rho} \der{p}{x_i}
   +
   \nu \der{}{x_j} \left( u_i \der{u_i}{x_j} \right)
   -
   \nu \der{u_i}{x_j} \der{u_i}{x_j}
   +
   u_i g_i,

where :math:`k \equiv \frac{1}{2} u_i u_i`.

