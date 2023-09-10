###################################
Derivation: internal energy balance
###################################

.. include:: tilde_note.rst

***************
Internal energy
***************

The internal energy of the liquid in a control volume leads

.. math::

   \int_{V} \rho c T dV \,\, \left[ M L^2 T^{-2} \right],

where :math:`c` is the specific heat capacity of the liquid.

The balance of this quantity is written as

.. math::

   \der{}{t} \int_{V} \rho c T dV
   =
   -
   \int_{V} \der{u_i \rho c T}{x_i} dV
   -
   \int_{\partial V} q_i n_i dS, \,\, \left[ M L^2 T^{-3} \right]

where :math:`q_i` is the heat flux given by

.. math::

   q_i
   =
   -
   k
   \der{T}{x_i}

using the `Fourier's law <https://en.wikipedia.org/wiki/Thermal_conduction#Fourier%27s_law>`_, where :math:`k` is the thermal conductivity.
Also, using the Gauss theorem, the diffusive term leads to

.. math::

   \int_{V} \der{}{x_i} \left( k \der{T}{x_i} \right) dV.

Assuming that the physical properties are constant, I obtain the internal energy balance:

.. math::

   \der{T}{t}
   +
   u_i \der{T}{x_i}
   =
   \kappa \der{}{x_i} \der{T}{x_i},

where :math:`\kappa` is the thermal diffusivity defined as :math:`k / \rho / c`.

*******************
Squared temperature
*******************

I consider to multiply the equation of the internal energy by :math:`T`.

* Temporal derivative

   .. math::

      T \der{T}{t}
      =
      \der{T^2}{t}
      -
      \der{T}{t} T,

   giving

   .. math::

      T \der{T}{t}
      =
      \frac{1}{2} \der{T^2}{t}.

* Advective terms

   .. math::

      T u_i \der{T}{x_i}
      =
      u_i \der{T^2}{x_i}
      -
      u_i \der{T}{x_i} T,

   giving

   .. math::

      T u_i \der{T}{x_i}
      =
      \frac{1}{2} u_i \der{T^2}{x_i}.

* Diffusive terms

   .. math::

      T \kappa \der{}{x_i} \der{T}{x_i}
      =
      \kappa \der{}{x_i} \left( T \der{T}{x_i} \right)
      -
      \kappa \der{T}{x_i} \der{T}{x_i}.

   .. note::

      The first term in the right-hand side yields

      .. math::

         \kappa \der{}{x_i} \left\{ \der{}{x_i} \left( T^2 \right) \right\}
         -
         \kappa \der{}{x_i} \left( \der{T}{x_i} T \right),

      and thus I notice

      .. math::

         \kappa \der{}{x_i} \left\{ \der{}{x_i} \left( \frac{1}{2} T^2 \right) \right\},

      which describes the diffusive effect of the squared temperature.

Thus, I have

.. math::

   \der{h}{t}
   +
   u_i \der{h}{x_i}
   =
   \kappa \der{}{x_i} \left( T \der{T}{x_i} \right)
   -
   \kappa \der{T}{x_i} \der{T}{x_i},

where :math:`h \equiv \frac{1}{2} T^2`.

