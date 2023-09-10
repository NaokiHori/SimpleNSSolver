Again I focus on the thermal energy equation; in particular I consider the equation of the **total** thermal energy :math:`\der{H}{t}`.
As discussed above, this quantity is increased by the thermal conduction on the walls:

.. math::

   & \int \frac{1}{\sqrt{Pr} \sqrt{Ra}} \der{}{\xi^j} \left( \mst{j}{j} T \der{T}{\xi^j} \right) d\gx d\gy \\
   \approx
   & \sum_i \sum_j \frac{1}{\sqrt{Pr} \sqrt{Ra}} \left\{
       \vat{\left( \mst{x}{x} T \diffe{T}{\gx} \right)}{\pip,\pjc}
     - \vat{\left( \mst{x}{x} T \diffe{T}{\gx} \right)}{\pim,\pjc}
     + \vat{\left( \mst{y}{y} T \diffe{T}{\gy} \right)}{\pic,\pjp}
     - \vat{\left( \mst{y}{y} T \diffe{T}{\gy} \right)}{\pic,\pjm}
   \right\} \\
   =
   & \sum_i \sum_j \frac{1}{\sqrt{Pr} \sqrt{Ra}} \left\{
       \vat{\left( \mst{x}{x} T \diffe{T}{\gx} \right)}{\pip,\pjc}
     - \vat{\left( \mst{x}{x} T \diffe{T}{\gx} \right)}{\pim,\pjc}
   \right\} \\
   =
   & \sum_j \frac{1}{\sqrt{Pr} \sqrt{Ra}} \left\{
       \vat{\left( \mst{x}{x} T \diffe{T}{\gx} \right)}{\text{right wall},\pjc}
     - \vat{\left( \mst{x}{x} T \diffe{T}{\gx} \right)}{\text{left  wall},\pjc}
   \right\}, \\
   =
   & - \sum_j \frac{1}{\sqrt{Pr} \sqrt{Ra}} \vat{\left( \frac{\Delta y_j}{\Delta x_{\text{left wall}}} \diffe{T}{\gx} \right)}{\text{left wall},\pjc},

where :math:`T_{\text{left wall}} \equiv 1` and :math:`T_{\text{right wall}} \equiv 0` are used.

On the other hand, :math:`H` is reduced by the dissipation:

.. math::

   & \int \frac{1}{\sqrt{Pr} \sqrt{Ra}} \der{}{\xi^j} \left( \mst{j}{j} T \der{T}{\xi^j} \right) d\gx d\gy \\
   \approx
   & \sum_i \sum_j \frac{1}{\sqrt{Pr} \sqrt{Ra}} \frac{1}{J_{\pic,\pjc}} \left\{
        \frac{1}{2} J_{\pip,\pjc} \left( \vat{\frac{\diffe{T}{\gx}}{\Delta x}}{\pip,\pjc} \right)^2
      + \frac{1}{2} J_{\pim,\pjc} \left( \vat{\frac{\diffe{T}{\gx}}{\Delta x}}{\pim,\pjc} \right)^2
      + \frac{1}{2} J_{\pic,\pjp} \left( \vat{\frac{\diffe{T}{\gy}}{\Delta y}}{\pic,\pjp} \right)^2
      + \frac{1}{2} J_{\pic,\pjm} \left( \vat{\frac{\diffe{T}{\gy}}{\Delta y}}{\pic,\pjm} \right)^2
   \right\}.

Again, I average this equation in the whole domain and in time.

The evolution of :math:`H` vanishes as usual since I am interested in the statistically-steady state.

The thermal conduction leads

.. math::

   - \frac{1}{l_x l_y} \sum_j \frac{1}{\sqrt{Pr} \sqrt{Ra}} \vat{\left( \frac{\Delta y_j}{\Delta x_{\text{left wall}}} \diffe{T}{\gx} \right)}{\text{left wall},\pjc} = Nu \,\, \left( \because l_x \equiv 1 \right),

while the dissipation leads

.. math::

   \frac{1}{l_x l_y} \sum_i \sum_j \frac{1}{\sqrt{Pr} \sqrt{Ra}} \left\{
        \frac{1}{2} J_{\pip,\pjc} \left( \vat{\frac{\diffe{T}{\gx}}{\Delta x}}{\pip,\pjc} \right)^2
      + \frac{1}{2} J_{\pim,\pjc} \left( \vat{\frac{\diffe{T}{\gx}}{\Delta x}}{\pim,\pjc} \right)^2
      + \frac{1}{2} J_{\pic,\pjp} \left( \vat{\frac{\diffe{T}{\gy}}{\Delta y}}{\pic,\pjp} \right)^2
      + \frac{1}{2} J_{\pic,\pjm} \left( \vat{\frac{\diffe{T}{\gy}}{\Delta y}}{\pic,\pjm} \right)^2
   \right\},

which is the thermal dissipation rate :math:`\ave{\epsilon_h}{x,y,t} = \ave{\der{T}{x_i} \der{T}{x_i}}{S,t}` in the continuous level.

Thus I notice

.. math::

   Nu = \sqrt{Pr} \sqrt{Ra} \ave{\epsilon_h}{x,y,t}.

