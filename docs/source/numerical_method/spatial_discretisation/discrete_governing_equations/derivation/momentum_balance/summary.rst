The following discrete equations are directly implemented in the code.

Momentum balance in the :math:`x` direction at :math:`\left( \xic, \xjc \right)`:

.. math::

   \der{\ux}{t}
   +
   \dintrpv{
      \dintrpa{\ux}{x}
      \dder{\ux}{x}
   }{x}
   +
   \dintrpa{
      \dintrpv{\uy}{x}
      \dder{\ux}{y}
   }{y}
   =
   -\dder{p}{x}
   +
   \frac{\sqrt{Pr}}{\sqrt{Ra}} \left\{
      \dder{}{x} \left( \dder{\ux}{x} \right)
      +
      \dder{}{y} \left( \dder{\ux}{y} \right)
   \right\}
   +
   \dintrpu{T}{x}.

Momentum balance in the :math:`y` direction at :math:`\left( \yic, \yjc \right)`:

.. math::

   \der{\uy}{t}
   +
   \dintrpv{
      \dintrpa{\ux}{y}
      \dder{\uy}{x}
   }{x}
   +
   \dintrpa{
      \dintrpa{\uy}{y}
      \dder{\uy}{y}
   }{y}
   =
   -\dder{p}{y}
   +
   \frac{\sqrt{Pr}}{\sqrt{Ra}} \left\{
      \dder{}{x} \left( \dder{\uy}{x} \right)
      +
      \dder{}{y} \left( \dder{\uy}{y} \right)
   \right\}.

.. seealso::

   :ref:`src/fluid/integrate/compute_rhs.c <fluid_compute_rhs>`.

.. note::

   I cannot define the thermal forcing term in the :math:`x` direction :math:`\dintrpu{T}{x}` here (recall that the *widetilde* symbol is a placeholder).
   To decide the correct interpolation, I need to take into account the equation of the internal energy and :ref:`the Nusselt number relations <nu_kinetic_energy_dissipation_discrete>`.

