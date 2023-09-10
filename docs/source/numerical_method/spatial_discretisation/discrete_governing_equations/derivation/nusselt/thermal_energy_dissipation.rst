
.. _nu_thermal_energy_dissipation_discrete:

######################################################
Nusselt number based on the thermal energy dissipation
######################################################

.. seealso::

   :ref:`Continuous counterpart <nu_thermal_energy_dissipation>`.

**********
Derivation
**********

I focus on :ref:`the discrete and volume-integrated equation of the squared temperature <eq_squared_temperature_discrete>`.

.. math::

   &
   \der{}{t} \left(
      \sum_{\pjc} \sum_{\pic} \frac{1}{2} T^2 \Delta x \Delta y
   \right) \\
   =
   &
   -
   \sum_{\pjc} \sum_{\pic}
   \frac{1}{\sqrt{Pr} \sqrt{Ra}}
   \left[
      \dintrpv{
         \left\{
            C
            \left( \dder{T}{x} \right)^2
         \right\}
      }{x}
      +
      \dintrpa{
         \left( \dder{T}{y} \right)^2
      }{y}
   \right]
   \Delta x \Delta y \\
   &
   +
   \sum_{\pjc}
   \frac{1}{\sqrt{Pr} \sqrt{Ra}}
   \left(
      \vat{
         T \dder{T}{x}
      }{x = 1}
      -
      \vat{
         T \dder{T}{x}
      }{x = 0}
   \right)
   \Delta y.

Since I am interested in the statistically-steady state, the left-hand-side term is zero and I have

.. math::

   &
   \sum_{\pjc} \sum_{\pic}
   \frac{1}{\sqrt{Pr} \sqrt{Ra}}
   \ave{
      \dintrpv{
         \left\{
            C
            \left( \dder{T}{x} \right)^2
         \right\}
      }{x}
      +
      \dintrpa{
         \left( \dder{T}{y} \right)^2
      }{y}
   }{t}
   \Delta x \Delta y \\
   &
   =
   \sum_{\pjc}
   \frac{1}{\sqrt{Pr} \sqrt{Ra}}
   \ave{
      \vat{
         T \dder{T}{x}
      }{x = 1}
      -
      \vat{
         T \dder{T}{x}
      }{x = 0}
   }{t}
   \Delta y.

The right-hand side yields

.. math::

   &
   \vat{T}{x = 1}
   \sum_{\pjc}
   \frac{1}{\sqrt{Pr} \sqrt{Ra}}
   \ave{
      \vat{
         \dder{T}{x}
      }{x = 1}
   }{t}
   \Delta y
   -
   \vat{T}{x = 0}
   \sum_{\pjc}
   \frac{1}{\sqrt{Pr} \sqrt{Ra}}
   \ave{
      \vat{
         T \dder{T}{x}
      }{x = 0}
   }{t}
   \Delta y \\
   =
   &
   -
   \vat{T}{x = 1}
   Nu J_{T,ref}^D
   +
   \vat{T}{x = 0}
   Nu J_{T,ref}^D \\
   =
   &
   Nu J_{T,ref}^D,

and thus

.. math::

   Nu
   =
   \frac{
      \sum_{\pjc} \sum_{\pic}
      \frac{1}{\sqrt{Pr} \sqrt{Ra}}
      \ave{
         \dintrpv{
            \left\{
               C
               \left( \dder{T}{x} \right)^2
            \right\}
         }{x}
         +
         \dintrpa{
            \left( \dder{T}{y} \right)^2
         }{y}
      }{t}
      \Delta x \Delta y
   }{
      J_{T,ref}^D
   }.

**************
Implementation
**************

I monitor the instantaneous value:

.. math::

   Nu \left( t \right)
   =
   \frac{
      \sum_{\pjc} \sum_{\pic}
      \frac{1}{\sqrt{Pr} \sqrt{Ra}}
      \left[
         \dintrpv{
            \left\{
               C
               \left( \dder{T}{x} \left( t \right) \right)^2
            \right\}
         }{x}
         +
         \dintrpa{
            \left( \dder{T}{y} \left( t \right) \right)^2
         }{y}
      \right]
      \Delta x \Delta y
   }{
      J_{T,ref}^D
   }.

All contributions in the numerator are added first, which is normalised by the reference value:

.. myliteralinclude:: /../../src/logging/nusselt/thermal_energy_dissipation.c
   :language: c
   :tag: compute thermal energy dissipation

* :math:`x` contribution

   .. myliteralinclude:: /../../src/logging/nusselt/thermal_energy_dissipation.c
      :language: c
      :tag: x contribution

* :math:`y` contribution

   .. myliteralinclude:: /../../src/logging/nusselt/thermal_energy_dissipation.c
      :language: c
      :tag: y contribution

* :math:`z` contribution

   .. myliteralinclude:: /../../src/logging/nusselt/thermal_energy_dissipation.c
      :language: c
      :tag: z contribution

