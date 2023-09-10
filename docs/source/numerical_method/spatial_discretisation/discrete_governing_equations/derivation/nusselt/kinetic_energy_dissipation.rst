
.. _nu_kinetic_energy_dissipation_discrete:

######################################################
Nusselt number based on the kinetic energy dissipation
######################################################

.. seealso::

   :ref:`Continuous counterpart <nu_kinetic_energy_dissipation>`.

**********
Derivation
**********

I focus on :ref:`the discrete and volume-integrated equation of the squared velocity <eq_squared_velocity_discrete>`.

.. math::

   &
   \der{}{t} \left(
      \sum_{\xjc} \sum_{\xic} \frac{1}{2} \ux^2 \Delta x \Delta y
      +
      \sum_{\yjc} \sum_{\yic} \frac{1}{2} \uy^2 \Delta x \Delta y
   \right) \\
   =
   &
   -
   \sum_{\xjc} \sum_{\xic}
   \frac{\sqrt{Pr}}{\sqrt{Ra}}
   \left[
      \dintrpv{
         \left( \dder{\ux}{x} \right)^2
      }{x}
      +
      \dintrpa{
         \left( \dder{\ux}{y} \right)^2
      }{y}
   \right]
   \Delta x \Delta y \\
   &
   -
   \sum_{\yjc} \sum_{\yic}
   \frac{\sqrt{Pr}}{\sqrt{Ra}}
   \left[
      \dintrpv{
         \left\{
            C
            \left( \dder{\uy}{x} \right)^2
         \right\}
      }{x}
      +
      \dintrpa{
         \left( \dder{\uy}{y} \right)^2
      }{y}
   \right]
   \Delta x \Delta y \\
   &
   +
   \sum_{\xjc} \sum_{\xic}
   \ux \dintrpu{T}{x}
   \Delta x \Delta y.

Since I am interested in the statistically-steady state, the left-hand-side term is zero and I have

.. math::

   &
   \sum_{\xjc} \sum_{\xic}
   \frac{\sqrt{Pr}}{\sqrt{Ra}}
   \ave{
      \dintrpv{
         \left( \dder{\ux}{x} \right)^2
      }{x}
      +
      \dintrpa{
         \left( \dder{\ux}{y} \right)^2
      }{y}
   }{t}
   \Delta x \Delta y \\
   +
   &
   \sum_{\yjc} \sum_{\yic}
   \frac{\sqrt{Pr}}{\sqrt{Ra}}
   \ave{
      \dintrpv{
         \left\{
            C
            \left( \dder{\uy}{x} \right)^2
         \right\}
      }{x}
      +
      \dintrpa{
         \left( \dder{\uy}{y} \right)^2
      }{y}
   }{t}
   \Delta x \Delta y \\
   =
   &
   \sum_{\xjc} \sum_{\xic}
   \ave{\ux \dintrpu{T}{x}}{t}
   \Delta x \Delta y.

The right-hand-side term is the kinetic energy injection which is discussed in :ref:`the previous page <nu_kinetic_energy_injection_discrete>`.
To equate them, finally I notice that the unknown interpolation

.. math::

   \dintrpu{T}{x}

should be

.. math::

   \dintrpa{T}{x}.

Also, from :ref:`the previous page <nu_kinetic_energy_injection_discrete>`, I have

.. math::

   Nu
   =
   \frac{
      \sum_{\pjc} \sum_{\xic}
      \ave{
         \ux
         \dintrpa{T}{x}
      }{t}
      \Delta x
      \Delta y
   }{
      J_{T,ref}^D
   }
   +
   1,

and by replacing the kinetic energy injection by the kinetic energy dissipation, I obtain

.. math::

   Nu
   & =
   \frac{
      \sum_{\xjc} \sum_{\xic}
      \frac{\sqrt{Pr}}{\sqrt{Ra}}
      \ave{
         \dintrpv{
            \left( \dder{\ux}{x} \right)^2
         }{x}
         +
         \dintrpa{
            \left( \dder{\ux}{y} \right)^2
         }{y}
      }{t}
      \Delta x \Delta y
   }{
      J_{T,ref}^D
   } \\
   & +
   \frac{
      \sum_{\yjc} \sum_{\yic}
      \frac{\sqrt{Pr}}{\sqrt{Ra}}
      \ave{
         \dintrpv{
            \left\{
               C
               \left( \dder{\uy}{x} \right)^2
            \right\}
         }{x}
         +
         \dintrpa{
            \left( \dder{\uy}{y} \right)^2
         }{y}
      }{t}
      \Delta x \Delta y
   }{
      J_{T,ref}^D
   } \\
   & +
   1,

which is the discrete analogue of the Nusselt number based on the kinetic energy dissipation.

**************
Implementation
**************

I monitor the instantaneous value:

.. math::

   Nu \left( t \right)
   & =
   \frac{
      \sum_{\xjc} \sum_{\xic}
      \frac{\sqrt{Pr}}{\sqrt{Ra}}
      \left[
         \dintrpv{
            \left( \dder{\ux}{x} \left( t \right) \right)^2
         }{x}
         +
         \dintrpa{
            \left( \dder{\ux}{y} \left( t \right) \right)^2
         }{y}
      \right]
      \Delta x \Delta y
   }{
      J_{T,ref}^D
   } \\
   & +
   \frac{
      \sum_{\yjc} \sum_{\yic}
      \frac{\sqrt{Pr}}{\sqrt{Ra}}
      \left[
         \dintrpv{
            \left\{
               C
               \left( \dder{\uy}{x} \left( t \right) \right)^2
            \right\}
         }{x}
         +
         \dintrpa{
            \left( \dder{\uy}{y} \left( t \right) \right)^2
         }{y}
      \right]
      \Delta x \Delta y
   }{
      J_{T,ref}^D
   } \\
   & +
   1.

The numerator (dissipation of the discrete kinetic energy) is summed (for each velocity component, for each direction), which is normalised by the reference value and the laminar contribution is added:

.. myliteralinclude:: /../../src/logging/nusselt/kinetic_energy_dissipation.c
   :language: c
   :tag: compute kinetic energy dissipation

* :math:`\ux-x` contribution

   .. myliteralinclude:: /../../src/logging/nusselt/kinetic_energy_dissipation.c
      :language: c
      :tag: ux-x contribution

* :math:`\ux-y` contribution

   .. myliteralinclude:: /../../src/logging/nusselt/kinetic_energy_dissipation.c
      :language: c
      :tag: ux-y contribution

* :math:`\ux-z` contribution

   .. myliteralinclude:: /../../src/logging/nusselt/kinetic_energy_dissipation.c
      :language: c
      :tag: ux-z contribution

* :math:`\uy-x` contribution

   .. myliteralinclude:: /../../src/logging/nusselt/kinetic_energy_dissipation.c
      :language: c
      :tag: uy-x contribution

* :math:`\uy-y` contribution

   .. myliteralinclude:: /../../src/logging/nusselt/kinetic_energy_dissipation.c
      :language: c
      :tag: uy-y contribution

* :math:`\uy-z` contribution

   .. myliteralinclude:: /../../src/logging/nusselt/kinetic_energy_dissipation.c
      :language: c
      :tag: uy-z contribution

* :math:`\uz-x` contribution

   .. myliteralinclude:: /../../src/logging/nusselt/kinetic_energy_dissipation.c
      :language: c
      :tag: uz-x contribution

* :math:`\uz-y` contribution

   .. myliteralinclude:: /../../src/logging/nusselt/kinetic_energy_dissipation.c
      :language: c
      :tag: uz-y contribution

* :math:`\uz-z` contribution

   .. myliteralinclude:: /../../src/logging/nusselt/kinetic_energy_dissipation.c
      :language: c
      :tag: uz-z contribution

