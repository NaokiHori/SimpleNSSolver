
.. _nu_heat_flux_discrete:

############################
Heat flux and Nusselt number
############################

.. seealso::

   :ref:`Continuous counterpart <nu_heat_flux>`.

*********
Heat flux
*********

I consider :ref:`the discrete equation of the internal energy <eq_temperature_discrete>`:

.. math::

   \der{T}{t}
   +
   \dder{
      \ux
      \dintrpa{T}{x}
   }{x}
   +
   \dder{
      \uy
      \dintrpa{T}{y}
   }{y}
   =
   \frac{1}{\sqrt{Pr} \sqrt{Ra}} \left\{
      \dder{}{x} \left( \dder{T}{x} \right)
      +
      \dder{}{y} \left( \dder{T}{y} \right)
   \right\},

which is defined at the cell center :math:`\left( \pic, \pjc \right)`, whose time average yields

.. math::

   \ave{
      \dder{
         \ux
         \dintrpa{T}{x}
      }{x}
   }{t}
   & +
   \ave{
      \dder{
         \uy
         \dintrpa{T}{y}
      }{y}
   }{t} \\
   & =
   \frac{1}{\sqrt{Pr} \sqrt{Ra}} \left\{
      \ave{
         \dder{}{x} \left( \dder{T}{x} \right)
      }{t}
      +
      \ave{
         \dder{}{y} \left( \dder{T}{y} \right)
      }{t}
   \right\}.

Now I consider to integrate this equation in the homogeneous direction (:math:`y`), giving

.. math::

   \sum_{\pjc}
   \ave{
      \dder{
         \ux
         \dintrpa{T}{x}
      }{x}
   }{t}
   \Delta y
   =
   \sum_{\pjc}
   \frac{1}{\sqrt{Pr} \sqrt{Ra}}
   \ave{
      \dder{}{x} \left( \dder{T}{x} \right)
   }{t}
   \Delta y

since the discretely-conservative terms in the :math:`y` direction

.. math::

   \dder{}{y} \left( \cdots \right)

disappear.

I further integrate this equation in the wall-normal direction from one place :math:`x = x_0` to the other :math:`x = x_1`:

.. math::

   \sum_{\pjc} \sum_{\pic}
   \ave{
      \dder{
         \ux
         \dintrpa{T}{x}
      }{x}
   }{t}
   \Delta x
   \Delta y
   =
   \sum_{\pjc} \sum_{\pic}
   \frac{1}{\sqrt{Pr} \sqrt{Ra}}
   \ave{
      \dder{}{x} \left( \dder{T}{x} \right)
   }{t}
   \Delta x
   \Delta y,

which is

.. math::

   \sum_{\pjc}
   \ave{
      \vat{
         \ux
         \dintrpa{T}{x}
      }{x_1}
      -
      \vat{
         \ux
         \dintrpa{T}{x}
      }{x_0}
   }{t}
   \Delta y
   =
   \sum_{\pjc}
   \frac{1}{\sqrt{Pr} \sqrt{Ra}} \left\{
      \vat{
         \ave{
            \dder{T}{x}
         }{t}
      }{x_1}
      -
      \vat{
         \ave{
            \dder{T}{x}
         }{t}
      }{x_0}
   \right\}
   \Delta y,

or equivalently

.. math::

   & \sum_{\pjc}
   \left[
      \ave{
         \ux
         \dintrpa{T}{x}
      }{t}
      -
      \frac{1}{\sqrt{Pr} \sqrt{Ra}}
      \ave{
         \dder{T}{x}
      }{t}
   \right]_{x_1}
   \Delta y \\
   & = \\
   & \sum_{\pjc}
   \left[
      \ave{
         \ux
         \dintrpa{T}{x}
      }{t}
      -
      \frac{1}{\sqrt{Pr} \sqrt{Ra}}
      \ave{
         \dder{T}{x}
      }{t}
   \right]_{x_0}
   \Delta y.

Notice that, since the integrand is defined at cell centers, :math:`x_0` and :math:`x_1` should be the :math:`x` cell faces.

Finally, I define the discrete time-averaged net heat flux:

.. _eq_heat_flux_discrete:

.. math::

   \vat{J_{T}^D}{\xic}
   \equiv
   \sum_{\pjc}
   \left[
      \ave{
         \ux
         \dintrpa{T}{x}
      }{t}
      -
      \frac{1}{\sqrt{Pr} \sqrt{Ra}}
      \ave{
         \dder{T}{x}
      }{t}
   \right]_{\xic}
   \Delta y,

which is constant for all :math:`x` faces.

.. note::

   In the continuous space, the net heat flux is constant at any :math:`x`.
   In the discrete space, on the other hand, the above relation holds only at each :math:`x` cell face, i.e. there is no guarantee that the above relation holds at cell centers.
   Also the same interpolating and the differential operators as the equation of the internal energy should be used to compute the above relation to rigorously satisfy the relation.

**************
Nusselt number
**************

To define the Nusselt number in the discrete space, I consider the discrete convective heat flux, which is

.. math::

   \vat{J_{T,ref}^D}{\xic}
   \equiv
   \sum_{\pjc}
   \left[
      -
      \frac{1}{\sqrt{Pr} \sqrt{Ra}}
      \ave{
         \dder{T_{ref}}{x}
      }{t}
   \right]_{\xic}
   \Delta y,

where :math:`T_{ref}` is the temperature field without the convective effect and simply given by solving the discrete second-order ordinary differential equation

.. math::

   \dder{}{x} \dder{\vat{T_{ref}}{\pic}}{x} = 0

under the Dirichlet boundary conditions:

.. math::

   \vat{T_{ref}}{x = 0} - \vat{T_{ref}}{x = 1} = 1.

.. note::

   I normalise the governing equations such that the temperature difference between the walls is unity.
   See :ref:`the governing equations <governing_equations>`.

giving

.. math::

   \vat{T_{ref}}{\pic}
   =
   C - \vat{x}{\pic},

or

.. math::

   \vat{\dder{T_{ref}}{x}}{\xic}
   =
   -1,

and thus

.. math::

   \vat{J_{T,ref}^D}{\xic}
   =
   \sum_{\pjc}
   \frac{1}{\sqrt{Pr} \sqrt{Ra}}
   \Delta y
   =
   const.

Now I can define the discrete Nusselt number as the ratio of the time-averaged discrete heat flux defined at the :math:`x` cell face:

.. math::

   \vat{J_{T}^D}{\xic}

to the discrete reference value without convection:

.. math::

   J_{T,ref}^D
   =
   \sum_{\pjc}
   \frac{1}{\sqrt{Pr} \sqrt{Ra}}
   \Delta y,

namely

.. _eq_nu_definition_discrete:

.. math::

   Nu
   \equiv
   \frac{
      \vat{J_{T}^D}{\xic}
   }{
      J_{T,ref}^D
   }
   =
   \frac{
      \sum_{\pjc}
      \left[
         \ave{
            \ux
            \dintrpa{T}{x}
         }{t}
         -
         \frac{1}{\sqrt{Pr} \sqrt{Ra}}
         \ave{
            \dder{T}{x}
         }{t}
      \right]_{\xic}
      \Delta y
   }{
      \sum_{\pjc}
      \frac{1}{\sqrt{Pr} \sqrt{Ra}}
      \Delta y
   }.

Note again that this relation is only valid where :math:`\ux` is defined.

To compute this quantity easily, I consider the above equation at the walls, giving

.. math::

   Nu_{wall}
   =
   \frac{
      J_{T,wall}^D
   }{
      J_{T,ref}^D
   }
   =
   \frac{
      \sum_{\pjc}
      \left[
         -
         \frac{1}{\sqrt{Pr} \sqrt{Ra}}
         \ave{
            \dder{T}{x}
         }{t}
      \right]_{wall}
      \Delta y
   }{
      \sum_{\pjc}
      \frac{1}{\sqrt{Pr} \sqrt{Ra}}
      \Delta y
   }.

**************
Implementation
**************

Numerically, I monitor the instantaneous value

.. math::

   Nu_{wall} \left( t \right)
   =
   \frac{
      \sum_{\pjc}
      \left[
         -
         \frac{1}{\sqrt{Pr} \sqrt{Ra}}
         \dder{T}{x} \left( t \right)
      \right]_{wall}
      \Delta y
   }{
      J_{T,ref}^D
   }

by assuming that the temporal and the spatial treatments are commutative, since I cannot perform the time average when the simulation is running.

.. myliteralinclude:: /../../src/logging/nusselt/heat_flux.c
   :language: c
   :tag: heat flux on the walls

