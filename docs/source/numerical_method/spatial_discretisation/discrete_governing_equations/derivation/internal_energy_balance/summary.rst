The following discrete equation is directly implemented in the code.

.. math::

   \der{T}{t}
   +
   \dintrpv{
      \ux
      \dder{T}{x}
   }{x}
   +
   \dintrpa{
      \uy
      \dder{T}{y}
   }{y}
   =
   \frac{1}{\sqrt{Pr} \sqrt{Ra}} \left\{
      \dder{}{x} \left( \dder{T}{x} \right)
      +
      \dder{}{y} \left( \dder{T}{y} \right)
   \right\},

which is defined at :math:`\left( \pic, \pjc \right)`.

.. seealso::

   :ref:`src/fluid/integrate/compute_rhs.c <fluid_compute_rhs>`.

