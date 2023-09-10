
.. _nu_kinetic_energy_dissipation:

######################################################
Nusselt number based on the kinetic energy dissipation
######################################################

I consider to integrate :ref:`the equation of the squared velocity <eq_squared_velocity>` in the whole volume.

* Temporal derivative

   .. math::

      \int_V \der{\ave{k}{t}}{t} dV

   is dropped because of the definition of the statistically-steady state.

* Advective and pressure-gradient contributions

   .. math::

      \int_V \ave{u_i \der{}{x_i} \left( - k - p \right)}{t} dV,

   which is

   .. math::

      \int_V \ave{\der{}{x_i} \left[ u_i \left( - k - p \right) \right]}{t} dV
      -
      \int_V \ave{\left( - k - p \right) \der{u_i}{x_i}}{t} dV,

   which also disappear because the first term is conservative and the second term includes the incompressibility.

* Diffusive contribution

   .. math::

      \int_V \frac{\sqrt{Pr}}{\sqrt{Ra}} \ave{\der{}{x_j} \left( u_i \der{u_i}{x_j} \right)}{t} dV
      -
      \int_V \frac{\sqrt{Pr}}{\sqrt{Ra}} \ave{\der{u_i}{x_j} \der{u_i}{x_j}}{t} dV,

   where the second term remains.
   Regarding the first term, the :math:`y` and :math:`z` derivatives are conservative and thus zero, whilst the :math:`x` derivative yields

   .. math::

      \int_V \frac{\sqrt{Pr}}{\sqrt{Ra}} \ave{\der{}{x} \left( u_i \der{u_i}{x} \right)}{t} dV
      & =
      \int_z \int_y \frac{\sqrt{Pr}}{\sqrt{Ra}} \ave{\left[ u_i \der{u_i}{x} \right]_{x = 1}}{t} dy dz \\
      & -
      \int_z \int_y \frac{\sqrt{Pr}}{\sqrt{Ra}} \ave{\left[ u_i \der{u_i}{x} \right]_{x = 0}}{t} dy dz.

   The integrands are all zero because :math:`u_i \equiv 0`, i.e.

   * :math:`\ux`: the walls are impermeable,

   * :math:`\uy`, :math:`\uz`: the walls are at rest.

   .. note::

      The second term gives the dissipation of the squared velocity (sink term).
      Since I can write the first term as

      .. math::

         \int_V \frac{1}{\sqrt{Pr} \sqrt{Ra}} \ave{\der{}{x_i} \left( \der{k}{x_i} \right)}{t} dV,

      I notice that this term tells the diffusive process of :math:`k`.
      In particular, this term is responsible for the transportation of :math:`k` between the walls and the fluid.

* Energy injection

   .. math::

      \int_V \ave{u_i T \delta_{ix}}{t} dV
      =
      \int_V \ave{\ux T}{t} dV.

In summary, I obtain

.. math::

   0
   =
   -
   \int_V \frac{\sqrt{Pr}}{\sqrt{Ra}} \ave{\der{u_i}{x_j} \der{u_i}{x_j}}{t} dV
   +
   \int_V \ave{\ux T}{t} dV.

As derived in :ref:`the previous page <nu_kinetic_energy_injection>`, I know

.. math::

   \int_V \ave{\ux T}{t} dV
   =
   J_{T,ref} \left( Nu - 1 \right).

Thus, I finally notice

.. math::

   \int_V \frac{\sqrt{Pr}}{\sqrt{Ra}} \ave{\der{u_i}{x_j} \der{u_i}{x_j}}{t} dV
   =
   J_{T,ref} \left( Nu - 1 \right),

or

.. math::

   Nu
   =
   \frac{1}{J_{T,ref}} \int_V \frac{\sqrt{Pr}}{\sqrt{Ra}} \ave{\der{u_i}{x_j} \der{u_i}{x_j}}{t} dV
   +
   1.

By introducing the local kinetic energy dissipation

.. math::

   \epsilon_k \left( x, y, z, t \right)
   \equiv
   \frac{\sqrt{Pr}}{\sqrt{Ra}} \der{u_i}{x_j} \der{u_i}{x_j},

I obtain the conclusive equation of this page:

.. math::

   Nu
   =
   \frac{1}{J_{T,ref}} \int_V \ave{\epsilon_k}{t} dV
   +
   1.

.. note::

   Since

   .. math::

      J_{T,ref}
      =
      \int_z \int_y \frac{1}{\sqrt{Pr} \sqrt{Ra}} dy dz
      =
      \int_V \frac{1}{\sqrt{Pr} \sqrt{Ra}} dV,

   I have

   .. math::

      Nu
      & =
      \frac{
         \int_V \ave{\epsilon_k}{t} dV
      }{
         \int_V \frac{1}{\sqrt{Pr} \sqrt{Ra}} dV
      }
      +
      1 \\
      & =
      \sqrt{Pr} \sqrt{Ra}
      \frac{
         \int_V \ave{\epsilon_k}{t} dV
      }{
         \int_V dV
      }
      +
      1 \\
      & =
      \sqrt{Pr} \sqrt{Ra} \ave{\epsilon_k}{V,t}
      +
      1,

   which is more generally used.

.. seealso::

   :ref:`Discrete counterpart <nu_kinetic_energy_dissipation_discrete>`.

