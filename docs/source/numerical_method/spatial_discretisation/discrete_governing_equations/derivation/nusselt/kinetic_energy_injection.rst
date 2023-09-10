
.. _nu_kinetic_energy_injection_discrete:

####################################################
Nusselt number based on the kinetic energy injection
####################################################

.. seealso::

   :ref:`Continuous counterpart <nu_kinetic_energy_injection>`.

**********
Derivation
**********

Recall that :ref:`the definition of the discrete Nusselt number <eq_nu_definition_discrete>` is

.. math::

   \vat{J_{T}^D}{\xic}
   =
   Nu
   \times
   J_{T,ref}^D.

Now I consider to integrate this equation from the bottom to the top walls.

The right-hand-side term yields

.. math::

   \sum_{\xic}
   Nu
   \times
   J_{T,ref}^D
   \Delta x
   & =
   Nu
   \times
   J_{T,ref}^D
   \sum_{\xic}
   \Delta x \\
   & =
   Nu
   \times
   J_{T,ref}^D \\
   & \left( \because l_x \equiv 1 \right),

while the left-hand side leads to

.. math::

   \sum_{\xic}
   J_{T}^D
   \Delta x
   & =
   \sum_{\xic}
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
   \right]
   \Delta y
   \Delta x \\
   & =
   \sum_{\pjc} \sum_{\xic}
   \ave{
      \ux
      \dintrpa{T}{x}
   }{t}
   \Delta x
   \Delta y \\
   & -
   \sum_{\pjc}
   \frac{1}{\sqrt{Pr} \sqrt{Ra}}
   \left(
      \ave{\vat{T}{x = 1}}{t}
      -
      \ave{\vat{T}{x = 0}}{t}
   \right)
   \Delta y.

The second term yields

.. math::

   +
   \sum_{\pjc}
   \frac{1}{\sqrt{Pr} \sqrt{Ra}}
   \Delta y

because of the boundary condition of the temperature field, and this is equal to :math:`J_{T,ref}^D`.

Thus I notice

.. math::

   Nu
   \times
   J_{T,ref}^D
   =
   \sum_{\pjc} \sum_{\xic}
   \ave{
      \ux
      \dintrpa{T}{x}
   }{t}
   \Delta x
   \Delta y
   +
   J_{T,ref}^D,

or equivalently


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

which is the discrete analogue of the Nusselt number based on the kinetic energy injection.

**************
Implementation
**************

I monitor the instantaneous value:

.. math::

   Nu \left( t \right)
   =
   \frac{
      \sum_{\pjc} \sum_{\xic}
      \ux \left( t \right)
      \dintrpa{T}{x} \left( t \right)
      \Delta x
      \Delta y
   }{
      J_{T,ref}^D
   }
   +
   1.

.. myliteralinclude:: /../../src/logging/nusselt/kinetic_energy_injection.c
   :language: c
   :tag: energy injection

