##########################################
:ref:`Momentum balance in y <eq_momentum>`
##########################################

Defined at :math:`\left( \yic, \yjc, \ykc \right)`:

.. math::

   \der{\uy}{t}
   =
   &
   -
   \dintrpv{
      \dintrpa{\ux}{y}
      \dder{\uy}{x}
   }{x}
   -
   \dintrpa{
      \dintrpa{\uy}{y}
      \dder{\uy}{y}
   }{y}
   -
   \dintrpa{
      \dintrpa{\uz}{y}
      \dder{\uy}{z}
   }{z} \\
   &
   -
   \dder{p}{y} \\
   &
   +
   \frac{\sqrt{Pr}}{\sqrt{Ra}} \left\{
      \dder{}{x} \left( \dder{\uy}{x} \right)
      +
      \dder{}{y} \left( \dder{\uy}{y} \right)
      +
      \dder{}{z} \left( \dder{\uy}{z} \right)
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

