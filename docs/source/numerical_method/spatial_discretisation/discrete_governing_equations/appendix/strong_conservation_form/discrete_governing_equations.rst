############################
Discrete governing equations
############################

The governing equations in the strong conservation form derived in the previous section are directly discretised based on the following notes.
For simplicity equations which are not treated numerically (kinetic energy, thermal quadratic quantity) are omitted.

.. note::

   * Since the new cooridnate system is designed to satisfy :math:`\Delta \xi \equiv \Delta \eta \equiv 1`, interpolations and differentiations are given as the arithmetic averages and the subtractions of two neighbouring values, respectively.

   * The Jacobian and the scaling factors should be interpolated as well as the variables.

   * The mesh skewness tensors lead to

      .. math::

         \mst{x}{x} = \frac{\Delta y}{\Delta x},
         \mst{y}{y} = \frac{\Delta x}{\Delta y}.

****************************
Incompressibility constraint
****************************

.. math::

   \frac{
       \vat{\ux}{\pip,\pjc}
     - \vat{\ux}{\pim,\pjc}
   }{\Delta x_{\pic}}
   + \frac{
       \vat{\uy}{\pic,\pjp}
     - \vat{\uy}{\pic,\pjm}
   }{\Delta y_{\pjc}}
   = 0,

which is evaluated at :math:`\left( \pic, \pjc \right)`.

****************
Momentum balance
****************

.. math::

   \der{\vat{\ux}{\xic,\xjc}}{t}
   =
   & - \frac{
       \vat{\dintrpa{\Delta y \ux}{\gx}}{\xip,\xjc} \vat{\dintrpa{\ux}{\gx}}{\xip,\xjc}
     - \vat{\dintrpa{\Delta y \ux}{\gx}}{\xim,\xjc} \vat{\dintrpa{\ux}{\gx}}{\xim,\xjc}
   }{\Delta x_{\xic} \Delta y_{\xjc}} \\
   & - \frac{
       \vat{\dintrpa{\Delta x \uy}{\gx}}{\xic,\xjp} \vat{\dintrpa{\ux}{\gy}}{\xic,\xjp}
     - \vat{\dintrpa{\Delta x \uy}{\gx}}{\xic,\xjm} \vat{\dintrpa{\ux}{\gy}}{\xic,\xjm}
   }{\Delta x_{\xic} \Delta y_{\xjc}} \\
   & -\frac{
       \vat{\Delta y p}{\xip,\xjc}
     - \vat{\Delta y p}{\xim,\xjc}
   }{\Delta x_{\xic} \Delta y_{\xjc}} \\
   & + \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{
       \vat{\frac{\Delta y}{\Delta x} \diffe{\ux}{\gx}}{\xip,\xjc}
     - \vat{\frac{\Delta y}{\Delta x} \diffe{\ux}{\gx}}{\xim,\xjc}
   }{\Delta x_{\xic} \Delta y_{\xjc}} \\
   & + \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{
       \vat{\frac{\Delta x}{\Delta y} \diffe{\ux}{\gy}}{\xic,\xjp}
     - \vat{\frac{\Delta x}{\Delta y} \diffe{\ux}{\gy}}{\xic,\xjm}
   }{\Delta x_{\xic} \Delta y_{\xjc}} \\
   & + \vat{\dintrpa{T}{\gx}}{\xic,\xjc},

which is evaluated at :math:`\left( \xic, \xjc \right)`.

.. math::

   \der{\vat{\uy}{\yic,\yjc}}{t}
   =
   & - \frac{
       \vat{\dintrpa{\Delta y \ux}{\gy}}{\yip,\yjc} \vat{\dintrpa{\uy}{\gx}}{\yip,\yjc}
     - \vat{\dintrpa{\Delta y \ux}{\gy}}{\yim,\yjc} \vat{\dintrpa{\uy}{\gx}}{\yim,\yjc}
   }{\Delta x_{\yic} \Delta y_{\yjc}} \\
   & - \frac{
       \vat{\dintrpa{\Delta x \uy}{\gy}}{\yic,\yjp} \vat{\dintrpa{\uy}{\gy}}{\yic,\yjp}
     - \vat{\dintrpa{\Delta x \uy}{\gy}}{\yic,\yjm} \vat{\dintrpa{\uy}{\gy}}{\yic,\yjm}
   }{\Delta x_{\yic} \Delta y_{\yjc}} \\
   & -\frac{
       \vat{\Delta x p}{\yic,\yjp}
     - \vat{\Delta x p}{\yic,\yjm}
   }{\Delta x_{\yic} \Delta y_{\yjc}} \\
   & + \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{
       \vat{\frac{\Delta y}{\Delta x} \diffe{\uy}{\gx}}{\yip,\yjc}
     - \vat{\frac{\Delta y}{\Delta x} \diffe{\uy}{\gx}}{\yim,\yjc}
   }{\Delta x_{\yic} \Delta y_{\yjc}} \\
   & + \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{
       \vat{\frac{\Delta x}{\Delta y} \diffe{\uy}{\gy}}{\yic,\yjp}
     - \vat{\frac{\Delta x}{\Delta y} \diffe{\uy}{\gy}}{\yic,\yjm}
   }{\Delta x_{\yic} \Delta y_{\yjc}},

which is evaluated at :math:`\left( \yic, \yjc \right)`.

***************
Internal energy
***************

.. math::

   \der{\vat{T}{\pic,\pjc}}{t}
   =
   & - \frac{
       \vat{\Delta y \ux}{\pip,\pjc} \vat{\dintrpa{T}{\gx}}{\pip,\pjc}
     - \vat{\Delta y \ux}{\pim,\pjc} \vat{\dintrpa{T}{\gx}}{\pim,\pjc}
   }{\Delta x_{\pic} \Delta y_{\pjc}} \\
   & - \frac{
       \vat{\Delta x \uy}{\pic,\pjp} \vat{\dintrpa{T}{\gy}}{\pic,\pjp}
     - \vat{\Delta x \uy}{\pic,\pjm} \vat{\dintrpa{T}{\gy}}{\pic,\pjm}
   }{\Delta x_{\pic} \Delta y_{\pjc}} \\
   & + \frac{1}{\sqrt{Pr} \sqrt{Ra}} \frac{
       \vat{\frac{\Delta y}{\Delta x} \diffe{T}{\gx}}{\pip,\pjc}
     - \vat{\frac{\Delta y}{\Delta x} \diffe{T}{\gx}}{\pim,\pjc}
   }{\Delta x_{\pic} \Delta y_{\pjc}} \\
   & + \frac{1}{\sqrt{Pr} \sqrt{Ra}} \frac{
       \vat{\frac{\Delta x}{\Delta y} \diffe{T}{\gy}}{\pic,\pjp}
     - \vat{\frac{\Delta x}{\Delta y} \diffe{T}{\gy}}{\pic,\pjm}
   }{\Delta x_{\pic} \Delta y_{\pjc}}.

which is evaluated at :math:`\left( \pic, \pjc \right)`.

