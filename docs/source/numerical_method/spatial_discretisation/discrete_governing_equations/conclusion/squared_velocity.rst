
.. _eq_squared_velocity_discrete:

#############################################################
:ref:`Equation of the squared velocity <eq_squared_velocity>`
#############################################################

Although the following relation is not implemented, it plays an important role in the conservation properties of the discretised equations.

Defined as the volume-integrated quantity:

.. math::

   &
   \der{}{t} \left(
      \sum_{\xkc} \sum_{\xjc} \sum_{\xic} \frac{1}{2} \ux^2 \Delta x \Delta y \Delta z
      +
      \sum_{\ykc} \sum_{\yjc} \sum_{\yic} \frac{1}{2} \uy^2 \Delta x \Delta y \Delta z
      +
      \sum_{\zkc} \sum_{\zjc} \sum_{\zic} \frac{1}{2} \uy^2 \Delta x \Delta y \Delta z
   \right) \\
   =
   &
   -
   \sum_{\xkc} \sum_{\xjc} \sum_{\xic}
   \frac{\sqrt{Pr}}{\sqrt{Ra}}
   \left[
      \dintrpv{
         \left( \dder{\ux}{x} \right)^2
      }{x}
      +
      \dintrpa{
         \left( \dder{\ux}{y} \right)^2
      }{y}
      +
      \dintrpa{
         \left( \dder{\ux}{z} \right)^2
      }{z}
   \right]
   \Delta x \Delta y \Delta z \\
   &
   -
   \sum_{\ykc} \sum_{\yjc} \sum_{\yic}
   \frac{\sqrt{Pr}}{\sqrt{Ra}}
   \left[
      \dintrpv{
         C
         \left( \dder{\uy}{x} \right)^2
      }{x}
      +
      \dintrpa{
         \left( \dder{\uy}{y} \right)^2
      }{y}
      +
      \dintrpa{
         \left( \dder{\uy}{z} \right)^2
      }{z}
   \right]
   \Delta x \Delta y \Delta z \\
   &
   -
   \sum_{\zkc} \sum_{\zjc} \sum_{\zic}
   \frac{\sqrt{Pr}}{\sqrt{Ra}}
   \left[
      \dintrpv{
         C
         \left( \dder{\uz}{x} \right)^2
      }{x}
      +
      \dintrpa{
         \left( \dder{\uz}{y} \right)^2
      }{y}
      +
      \dintrpa{
         \left( \dder{\uz}{z} \right)^2
      }{z}
   \right]
   \Delta x \Delta y \Delta z \\
   &
   +
   \sum_{\xkc} \sum_{\xjc} \sum_{\xic}
   \ux \dintrpa{T}{x}
   \Delta x \Delta y \Delta z,

where :math:`C` is a constant to take into account the no-slip boundary effect.

.. seealso::

   * :ref:`logging/nusselt/kinetic_energy_dissipation.c <nu_kinetic_energy_dissipation_discrete>`

