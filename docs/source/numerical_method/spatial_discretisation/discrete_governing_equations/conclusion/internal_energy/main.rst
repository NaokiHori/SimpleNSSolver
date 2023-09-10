.. _eq_temperature_discrete:

###############################################
:ref:`Internal energy balance <eq_temperature>`
###############################################

Defined at :math:`\left( \pic, \pjc, \pkc \right)`:

.. math::

   \der{T}{t}
   =
   &
   -
   \dintrpv{
      \ux
      \dder{T}{x}
   }{x}
   -
   \dintrpa{
      \uy
      \dder{T}{y}
   }{y}
   -
   \dintrpa{
      \uz
      \dder{T}{z}
   }{z} \\
   &
   +
   \frac{1}{\sqrt{Pr} \sqrt{Ra}} \left\{
      \dder{}{x} \left( \dder{T}{x} \right)
      +
      \dder{}{y} \left( \dder{T}{y} \right)
      +
      \dder{}{z} \left( \dder{T}{z} \right)
   \right\}.

.. toctree::
   :maxdepth: 1
   :caption: Implmentation

   adv_x
   adv_y
   adv_z
   dif_x
   dif_y
   dif_z

.. seealso::

   * :ref:`fluid_compute_rhs <fluid_compute_rhs>`

