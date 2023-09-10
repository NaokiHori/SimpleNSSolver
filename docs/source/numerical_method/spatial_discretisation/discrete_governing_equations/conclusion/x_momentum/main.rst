##########################################
:ref:`Momentum balance in x <eq_momentum>`
##########################################

Defined at :math:`\left( \xic, \xjc, \xkc \right)`:

.. math::

   \der{\ux}{t}
   =
   &
   -
   \dintrpv{
      \dintrpa{\ux}{x}
      \dder{\ux}{x}
   }{x}
   -
   \dintrpa{
      \dintrpv{\uy}{x}
      \dder{\ux}{y}
   }{y}
   -
   \dintrpa{
      \dintrpv{\uz}{x}
      \dder{\ux}{z}
   }{z} \\
   &
   -
   \dder{p}{x} \\
   &
   +
   \frac{\sqrt{Pr}}{\sqrt{Ra}} \left\{
      \dder{}{x} \left( \dder{\ux}{x} \right)
      +
      \dder{}{y} \left( \dder{\ux}{y} \right)
      +
      \dder{}{z} \left( \dder{\ux}{z} \right)
   \right\} \\
   &
   +
   \dintrpa{T}{x}.

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

