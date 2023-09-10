
.. _nu_thermal_energy_dissipation:

######################################################
Nusselt number based on the thermal energy dissipation
######################################################

I consider to integrate :ref:`the equation of the squared temperature <eq_squared_temperature>` in the whole volume.

* Temporal derivative

   .. math::

      \int_V \der{\ave{h}{t}}{t} dV

   is dropped because of the definition of the statistically-steady state.

* Advective contribution

   .. math::

      \int_V \ave{u_i \der{h}{x_i}}{t} dV,

   which is

   .. math::

      \int_V \ave{\der{}{x_i} \left( - u_i - h \right)}{t} dV
      -
      \int_V \ave{\left( - h \right) \der{u_i}{x_i}}{t} dV,

   and also disappears because the first term is conservative and the second term includes the incompressibility.

* Diffusive contribution

   .. math::

      \int_V \frac{1}{\sqrt{Pr} \sqrt{Ra}} \ave{\der{}{x_i} \left( T \der{T}{x_i} \right)}{t} dV
      -
      \int_V \frac{1}{\sqrt{Pr} \sqrt{Ra}} \ave{\der{T}{x_i} \der{T}{x_i}}{t} dV,

   where the second term remains as it is.
   Regarding the first term, the :math:`y` and :math:`z` derivatives are conservative and thus zero, whilst the :math:`x` component yields

   .. math::

      &
      \int_z \int_y \frac{1}{\sqrt{Pr} \sqrt{Ra}} \left[ \ave{T \der{T}{x}}{t} \right]_{x = 1} dy dz
      -
      \int_z \int_y \frac{1}{\sqrt{Pr} \sqrt{Ra}} \left[ \ave{T \der{T}{x}}{t} \right]_{x = 0} dy dz \\
      & =
      \vat{\ave{T}{t}}{x = 1} \int_z \int_y \frac{1}{\sqrt{Pr} \sqrt{Ra}} \left[ \ave{\der{T}{x}}{t} \right]_{x = 1} dy dz \\
      & -
      \vat{\ave{T}{t}}{x = 0} \int_z \int_y \frac{1}{\sqrt{Pr} \sqrt{Ra}} \left[ \ave{\der{T}{x}}{t} \right]_{x = 0} dy dz,

   where the integrands are :ref:`the heat fluxes <nu_heat_flux>` evaluated at :math:`x = 1` and at :math:`x = 0`, respectively, and thus

   .. math::

      -
      \vat{\ave{T}{t}}{x = 1} J_{T} \left( x = 1 \right)
      +
      \vat{\ave{T}{t}}{x = 0} J_{T} \left( x = 0 \right)
      =
      \left[
         -
         \vat{\ave{T}{t}}{x = 1}
         +
         \vat{\ave{T}{t}}{x = 0}
      \right]
      J_{T}.

   Since I know

   .. math::

      \vat{\ave{T}{t}}{x = 0}
      -
      \vat{\ave{T}{t}}{x = 1}
      \equiv
      1

   because of the boundary conditions and

   .. math::

      J_{T} \left( x \right)
      =
      J_{T}
      =
      Nu
      J_{T,ref}

   because of :ref:`the defnition of the Nusselt number <eq_nu_definition>`, I finally find that the diffusive term yields

   .. math::

      Nu J_{T,ref}
      -
      \int_V \frac{1}{\sqrt{Pr} \sqrt{Ra}} \ave{\der{T}{x_i} \der{T}{x_i}}{t} dV.

   .. note::

      The second term gives the dissipation of the squared temperature (sink term).
      Since I can write the first term as

      .. math::

         \int_V \frac{1}{\sqrt{Pr} \sqrt{Ra}} \ave{\der{}{x_i} \left( \der{h}{x_i} \right)}{t} dV,

      I notice that this term tells the diffusive process of :math:`h`.
      In particular, this term is responsible for the transportation of :math:`h` between the walls and the fluid.

In summary, I obtain

.. math::

   0
   =
   Nu J_{T,ref}
   -
   \int_V \frac{1}{\sqrt{Pr} \sqrt{Ra}} \ave{\der{T}{x_i} \der{T}{x_i}}{t} dV,

or

.. math::

   Nu
   =
   \frac{1}{J_{T,ref}}
   \int_V \frac{1}{\sqrt{Pr} \sqrt{Ra}} \ave{\der{T}{x_i} \der{T}{x_i}}{t} dV.

By introducing the local thermal energy dissipation

.. math::

   \epsilon_h \left( x, y, z, t \right)
   \equiv
   \frac{1}{\sqrt{Pr} \sqrt{Ra}} \der{T}{x_i} \der{T}{x_i},

I obtain the conclusive equation of this page:

.. math::

   Nu
   =
   \frac{1}{J_{T,ref}} \int_V \ave{\epsilon_h}{t} dV.

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
         \int_V \epsilon_h dV
      }{
         \int_V \frac{1}{\sqrt{Pr} \sqrt{Ra}} dV
      } \\
      & =
      \sqrt{Pr} \sqrt{Ra}
      \frac{
         \int_V \ave{\epsilon_h}{t} dV
      }{
         \int_V dV
      } \\
      & =
      \sqrt{Pr} \sqrt{Ra} \ave{\epsilon_h}{V,t},

   which is more generally used.

.. seealso::

   :ref:`Discrete counterpart <nu_thermal_energy_dissipation_discrete>`.

