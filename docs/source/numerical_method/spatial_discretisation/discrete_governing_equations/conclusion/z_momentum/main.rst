##########################################
:ref:`Momentum balance in z <eq_momentum>`
##########################################

Defined at :math:`\left( \zic, \zjc, \zkc \right)`:

.. math::

   \der{\uz}{t}
   =
   &
   -
   \dintrpv{
      \dintrpa{\ux}{z}
      \dder{\uz}{x}
   }{x}
   -
   \dintrpa{
      \dintrpa{\uy}{z}
      \dder{\uz}{y}
   }{y}
   -
   \dintrpa{
      \dintrpa{\uz}{z}
      \dder{\uz}{z}
   }{z} \\
   &
   -
   \dder{p}{z} \\
   &
   +
   \frac{\sqrt{Pr}}{\sqrt{Ra}} \left\{
      \dder{}{x} \left( \dder{\uz}{x} \right)
      +
      \dder{}{y} \left( \dder{\uz}{y} \right)
      +
      \dder{}{z} \left( \dder{\uz}{z} \right)
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

