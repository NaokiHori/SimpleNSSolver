
.. _eq_squared_temperature_discrete:

###################################################################
:ref:`Equation of the squared temperature <eq_squared_temperature>`
###################################################################

Although the following relation is not implemented, it plays an important role in the conservation properties of the discretised equations.

Defined as the volume-integrated quantity:

.. math::

   &
   \der{}{t} \left(
      \sum_{\pkc} \sum_{\pjc} \sum_{\pic} \frac{1}{2} T^2 \Delta x \Delta y \Delta z
   \right) \\
   =
   &
   -
   \sum_{\pkc} \sum_{\pjc} \sum_{\pic}
   \frac{1}{\sqrt{Pr} \sqrt{Ra}}
   \left[
      \dintrpv{
         C
         \left( \dder{T}{x} \right)^2
      }{x}
      +
      \dintrpa{
         \left( \dder{T}{y} \right)^2
      }{y}
      +
      \dintrpa{
         \left( \dder{T}{z} \right)^2
      }{z}
   \right]
   \Delta x \Delta y \Delta z \\
   &
   +
   \sum_{\pkc} \sum_{\pjc}
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
   \Delta y \Delta z,

where :math:`C` is a constant to take into account the no-slip boundary effect.

.. seealso::

   * :ref:`logging/nusselt/thermal_energy_dissipation.c <nu_thermal_energy_dissipation_discrete>`

