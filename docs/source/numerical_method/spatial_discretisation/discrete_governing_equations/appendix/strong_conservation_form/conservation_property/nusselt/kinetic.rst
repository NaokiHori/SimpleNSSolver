Next, I move on to the kinetic energy equation to find the other relation.
As discussed above, I have the equation of the total kinetic energy :math:`\der{K}{t}`, which is increased by the energy injection by the buoyancy force:

.. math::

   \sum_{\forall \ux \text{positions}, \left( \pip, \pjc \right)} J_{\pip,\pjc} \vat{\ux}{\pip,\pjc} \vat{\dintrpa{T}{\gx}}{\pip,\pjc},

while reduced by the viscous dissipation coming from the momentum equation in :math:`x` direction:

.. math::

   \sum_{\forall \ux \text{positions}, \left( \xic, \xjc \right)} \frac{\sqrt{Pr}}{\sqrt{Ra}} \left\{
        \frac{1}{2} J_{\xip,\xjc} \left( \vat{\frac{\diffe{\ux}{\gx}}{\Delta x}}{\xip,\xjc} \right)^2
      + \frac{1}{2} J_{\xim,\xjc} \left( \vat{\frac{\diffe{\ux}{\gx}}{\Delta x}}{\xim,\xjc} \right)^2
      + \frac{1}{2} J_{\xic,\xjp} \left( \vat{\frac{\diffe{\ux}{\gy}}{\Delta y}}{\xic,\xjp} \right)^2
      + \frac{1}{2} J_{\xic,\xjm} \left( \vat{\frac{\diffe{\ux}{\gy}}{\Delta y}}{\xic,\xjm} \right)^2
   \right\}

and the one in :math:`y` direction:

.. math::

   \sum_{\forall \uy \text{positions}, \left( \yic, \yjc \right)} \frac{\sqrt{Pr}}{\sqrt{Ra}} \left\{
        \frac{1}{2} J_{\yip,\yjc} \left( \vat{\frac{\diffe{\uy}{\gx}}{\Delta x}}{\yip,\yjc} \right)^2
      + \frac{1}{2} J_{\yim,\yjc} \left( \vat{\frac{\diffe{\uy}{\gx}}{\Delta x}}{\yim,\yjc} \right)^2
      + \frac{1}{2} J_{\yic,\yjp} \left( \vat{\frac{\diffe{\uy}{\gy}}{\Delta y}}{\yic,\yjp} \right)^2
      + \frac{1}{2} J_{\yic,\yjm} \left( \vat{\frac{\diffe{\uy}{\gy}}{\Delta y}}{\yic,\yjm} \right)^2
   \right\}.

Again, I consider to average this equation in the whole domain and in time.

The evolution of :math:`K` vanishes as usual since I am interested in the statistically-steady state.
The buoyancy forcing leads

.. math::

   \frac{1}{l_x l_y} \sum_{\forall \ux \text{positions}, \left( \pip, \pjc \right)} J_{\pip,\pjc} \vat{\ux}{\pip,\pjc} \vat{\dintrpa{T}{\gx}}{\pip,\pjc}.

The total viscous dissipation is given as the sum of

.. math::

   \frac{1}{l_x l_y} \sum_{\forall \ux \text{positions}, \left( \xic, \xjc \right)} \frac{\sqrt{Pr}}{\sqrt{Ra}} \left\{
        \frac{1}{2} J_{\xip,\xjc} \left( \vat{\frac{\diffe{\ux}{\gx}}{\Delta x}}{\xip,\xjc} \right)^2
      + \frac{1}{2} J_{\xim,\xjc} \left( \vat{\frac{\diffe{\ux}{\gx}}{\Delta x}}{\xim,\xjc} \right)^2
      + \frac{1}{2} J_{\xic,\xjp} \left( \vat{\frac{\diffe{\ux}{\gy}}{\Delta y}}{\xic,\xjp} \right)^2
      + \frac{1}{2} J_{\xic,\xjm} \left( \vat{\frac{\diffe{\ux}{\gy}}{\Delta y}}{\xic,\xjm} \right)^2
   \right\}

and

.. math::

   \frac{1}{l_x l_y} \sum_{\forall \uy \text{positions}, \left( \yic, \yjc \right)} \frac{\sqrt{Pr}}{\sqrt{Ra}} \left\{
        \frac{1}{2} J_{\yip,\yjc} \left( \vat{\frac{\diffe{\uy}{\gx}}{\Delta x}}{\yip,\yjc} \right)^2
      + \frac{1}{2} J_{\yim,\yjc} \left( \vat{\frac{\diffe{\uy}{\gx}}{\Delta x}}{\yim,\yjc} \right)^2
      + \frac{1}{2} J_{\yic,\yjp} \left( \vat{\frac{\diffe{\uy}{\gy}}{\Delta y}}{\yic,\yjp} \right)^2
      + \frac{1}{2} J_{\yic,\yjm} \left( \vat{\frac{\diffe{\uy}{\gy}}{\Delta y}}{\yic,\yjm} \right)^2
   \right\}.

Note that this summation is nothing else but the (discretised) kinetic dissipation rate :math:`\ave{\epsilon_k}{x,y,t}`, i.e. :math:`\left\langle s_{ij} s_{ij} \right\rangle_{S,t}` in the continuous level.

Finally I find

.. math::

   \frac{1}{l_x l_y} \sum_{\forall \ux \text{positions}, \left( \pip, \pjc \right)} J_{\pip,\pjc} \vat{\ux}{\pip,\pjc} \vat{\dintrpa{T}{\gx}}{\pip,\pjc}
   =
   \ave{\epsilon_k}{x,y,t},

and, by using the relation I obtained in the previous Nusselt definition:

.. math::

   Nu = \frac{1}{l_x l_y} \sum_{\forall \ux \text{positions}, \left( \pip, \pjc \right)} \left(
      \sqrt{Pr} \sqrt{Ra} \, J_{\pip,\pjc} \vat{\ux}{\pip,\pjc} \vat{\dintrpa{T}{\gx}}{\pip,\pjc}
   \right)
   + 1,

I obtain

.. math::

   Nu = \sqrt{Pr} \sqrt{Ra} \ave{\epsilon_k}{x,y,t} + 1.

