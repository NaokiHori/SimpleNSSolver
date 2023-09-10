
.. _nu_heat_flux:

############################
Heat flux and Nusselt number
############################

To quantitatively discuss the enhancement and to define :math:`Nu` eventually, I need an equation describing the heat flux between two walls (how much the thermal energy is transported).

*********
Heat flux
*********

To begin with, I look at the equation of the (non-dimensional) internal energy (temperature):

.. math::

   \der{T}{t}
   +
   u_i \der{T}{x_i}
   =
   \frac{1}{\sqrt{Pr} \sqrt{Ra}} \der{}{x_i} \der{T}{x_i},

and integrate this equation in the homogeneous directions (:math:`y` and :math:`z`), since I am now interested in the total heat flux.

Since I am interested in systems which are in statistically-steady states, I first consider to average the equation in time:

.. math::

   \ave{q}{t}
   \equiv
   \frac{1}{t_{max} - t_{min}}
   \int_{\tau = t_{min}}^{\tau = t_{max}} q d \tau,

giving

.. math::

   \ave{u_i \der{T}{x_i}}{t}
   =
   \frac{1}{\sqrt{Pr} \sqrt{Ra}} \ave{\der{}{x_i} \der{T}{x_i}}{t}.

.. note::

   Hereafter I assume that operations in time and in space are commutative, which is in general not accepted but hold in the current context.

By using the incompressibility constraint, the advective terms lead to

.. math::

   \der{u_i T}{x_i}
   =
   \der{\ux T}{x}
   +
   \der{\uy T}{y}
   +
   \der{\uz T}{z}.

Since they are conservative, the last two terms vanish because of the periodic boundary condition, e.g.

.. math::

   \int_{y = 0}^{y = l_y} \der{\uy T}{y} dy
   =
   \left[ \uy T \right]_{y = l_y}
   -
   \left[ \uy T \right]_{y =   0},

which is zero because of the periodic boundary condition:

.. math::

   \vat{q}{y = l_y}
   =
   \vat{q}{y = 0}.

The second-order derivatives in the :math:`y` and :math:`z` directions, which reflect the diffusive effect, are also zero for the same reason, and thus I am left with

.. math::

   \int_{z} \int_{y}
   \left(
      \ave{\der{\ux T}{x}}{t}
      -
      \frac{1}{\sqrt{Pr} \sqrt{Ra}} \ave{\der{}{x} \der{T}{x}}{t}
   \right)
   dy dz
   =
   0.

Then, I further integrate this equation in the wall-normal direction from one place :math:`x = x_0` to the other :math:`x = x_1`:

.. math::

   \int_{x = x_0}^{x = x_1} \int_{z} \int_{y}
   \left(
      \ave{\der{\ux T}{x}}{t}
      -
      \frac{1}{\sqrt{Pr} \sqrt{Ra}} \ave{\der{}{x} \der{T}{x}}{t}
   \right)
   dy dz dx
   =
   0,

or equivalently

.. math::

   \int_{z} \int_{y} \left[ \ave{\ux T}{t} - \frac{1}{\sqrt{Pr} \sqrt{Ra}} \ave{\der{T}{x}}{t} \right]_{x = x_0} dy dz
   =
   \int_{z} \int_{y} \left[ \ave{\ux T}{t} - \frac{1}{\sqrt{Pr} \sqrt{Ra}} \ave{\der{T}{x}}{t} \right]_{x = x_1} dy dz,

where :math:`x_0` and :math:`x_1` are arbitrary as long as :math:`\in \left[ 0, l_x \equiv 1 \right]`.
This gives the definition of the time-averaged net heat flux :math:`J_{T} \left( x \right)`:

.. _eq_heat_flux:

.. math::

   J_T \left( x \right)
   \equiv
   \int_{z} \int_{y} \left[
      \ave{\ux T}{t}
      -
      \frac{1}{\sqrt{Pr} \sqrt{Ra}} \ave{\der{T}{x}}{t}
   \right] dy dz
   =
   const.

Two important takeaways are as follows.

* The net heat flux is consisted of the advective effect (the first term) and the diffusive effect (the second term).

* The net heat flux in the wall-normal direction (heat flux going through an arbitrary :math:`y-z` plane inside the domain) is constant.

.. note::

   The integrand is the time-averaged and non-dimensional `differential form of Fourier's law <https://en.wikipedia.org/wiki/Thermal_conduction#Differential_form>`_ if the advective term is absent:

   .. math::

      - \frac{1}{\sqrt{Pr} \sqrt{Ra}} \ave{\der{T}{x}}{t}.

**************
Nusselt number
**************

Recall that I am interested in quantifying the heat transfer enhancement by the convective effect.
To do so, I define :math:`Nu` as

.. math::

   Nu
   \equiv
   \frac{J_{T}}{J_{T,ref}},

where the denominator is the *reference* net heat flux, which is the heat flux if the convective effect is missing.

If the convection is absent (:math:`u_i \equiv 0`), the time-averaged equation of the internal energy is simply

.. math::

   \der{}{x} \der{\ave{T}{t}}{x}
   =
   0,

with the boundary conditions

.. math::

   \vat{T}{x = 0}
   -
   \vat{T}{x = 1}
   =
   1,

which have an analytical solution:

.. math::

   T \left( x \right)
   =
   const. - x
   \,\, \rightarrow \,\,
   \der{T}{x}
   =
   -1,

and thus

.. math::

   J_{T,ref}
   =
   \int_{z} \int_{y} \left( - \frac{1}{\sqrt{Pr} \sqrt{Ra}} \times \left( -1 \right) \right) dy dz
   =
   \int_{z} \int_{y} \frac{1}{\sqrt{Pr} \sqrt{Ra}} dy dz.

In summary, I have

.. _eq_nu_definition:

.. math::

   Nu
   \equiv
   \frac{J_{T} \left( x \right)}{J_{T,ref}}
   =
   \frac{
      \int_{z} \int_{y} \left(
         \ave{\ux T}{t}
         -
         \frac{1}{\sqrt{Pr} \sqrt{Ra}} \ave{\der{T}{x}}{t}
      \right) dy dz
   }{
      \int_{z} \int_{y} \frac{1}{\sqrt{Pr} \sqrt{Ra}} dy dz
   },

which I regard as the **definition** of the Nusselt number in this project.

Note that the numerator :math:`J_{T} \left( x \right)` is constant at any :math:`x` and thus can be evaluated at any :math:`x`.
However, the most intuitive way would be to compute on the walls, i.e.

.. math::

   J_{T} \left( x = 0 \right)
   =
   \int_{z} \int_{y} \left( - \frac{1}{\sqrt{Pr} \sqrt{Ra}} \vat{\ave{\der{T}{x}}{t}}{x = 0} \right) dy dz

at :math:`x = 0`, or

.. math::

   J_{T} \left( x = 1 \right)
   =
   \int_{z} \int_{y} \left( - \frac{1}{\sqrt{Pr} \sqrt{Ra}} \vat{\ave{\der{T}{x}}{t}}{x = 1} \right) dy dz

at :math:`x = l_x \equiv 1`, where the advective term is dropped because of the impermeable condition :math:`\vat{\ux}{x = 0} = \vat{\ux}{x = 1} \equiv 0`.

Finally I obtain the following famous relation:

.. math::

   Nu
   =
   \frac{J_{T} \left( x = wall \right)}{J_{T,ref}}
   =
   \frac{
      \int_{z} \int_{y} \left( - \frac{1}{\sqrt{Pr} \sqrt{Ra}} \vat{\ave{\der{T}{x}}{t}}{wall} \right) dy dz
   }{
      \int_{z} \int_{y} \frac{1}{\sqrt{Pr} \sqrt{Ra}} dy dz
   }.

.. seealso::

   :ref:`Discrete counterpart <nu_heat_flux_discrete>`.

