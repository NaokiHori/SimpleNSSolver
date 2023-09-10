###################################
Discrete thermal quadratic quantity
###################################

The equation of the squared temperature leads to

.. math::

   \vat{T}{\pic,\pjc} \der{\vat{T}{\pic,\pjc}}{t}
   =
   & - \vat{T}{\pic,\pjc} \frac{
       \vat{\Delta y \ux}{\pip,\pjc} \vat{\dintrpa{T}{\gx}}{\pip,\pjc}
     - \vat{\Delta y \ux}{\pim,\pjc} \vat{\dintrpa{T}{\gx}}{\pim,\pjc}
   }{\Delta x_{\pic} \Delta y_{\pjc}} \\
   & - \vat{T}{\pic,\pjc} \frac{
       \vat{\Delta x \uy}{\pic,\pjp} \vat{\dintrpa{T}{\gy}}{\pic,\pjp}
     - \vat{\Delta x \uy}{\pic,\pjm} \vat{\dintrpa{T}{\gy}}{\pic,\pjm}
   }{\Delta x_{\pic} \Delta y_{\pjc}} \\
   & + \vat{T}{\pic,\pjc} \frac{1}{\sqrt{Pr} \sqrt{Ra}} \frac{
       \vat{\frac{\Delta y}{\Delta x} \diffe{T}{\gx}}{\pip,\pjc}
     - \vat{\frac{\Delta y}{\Delta x} \diffe{T}{\gx}}{\pim,\pjc}
   }{\Delta x_{\pic} \Delta y_{\pjc}} \\
   & + \vat{T}{\pic,\pjc} \frac{1}{\sqrt{Pr} \sqrt{Ra}} \frac{
       \vat{\frac{\Delta x}{\Delta y} \diffe{T}{\gy}}{\pic,\pjp}
     - \vat{\frac{\Delta x}{\Delta y} \diffe{T}{\gy}}{\pic,\pjm}
   }{\Delta x_{\pic} \Delta y_{\pjc}}.

***************
Advective terms
***************

The advection terms of the squared temperature leads to

.. math::

   & - \frac{
       \vat{\left( \Delta y \ux \right)}{\pip,\pjc} \frac{\vat{T}{\pic ,\pjc } \vat{T}{\pipp,\pjc }}{2}
     - \vat{\left( \Delta y \ux \right)}{\pim,\pjc} \frac{\vat{T}{\pimm,\pjc } \vat{T}{\pic ,\pjc }}{2}
   }{\Delta x_{\pic} \Delta y_{\pjc}} \\
   & - \frac{
       \vat{\left( \Delta x \uy \right)}{\pic,\pjp} \frac{\vat{T}{\pic ,\pjc } \vat{T}{\pic ,\pjpp}}{2}
     - \vat{\left( \Delta x \uy \right)}{\pic,\pjm} \frac{\vat{T}{\pic ,\pjmm} \vat{T}{\pic ,\pjc }}{2}
   }{\Delta x_{\pic} \Delta y_{\pjc}} \\
   & - \frac{\vat{T^2}{\pic,\pjc}}{2} \left(
      \frac{
          \vat{\Delta y \ux}{\pip,\pjc}
        - \vat{\Delta y \ux}{\pim,\pjc}
      }{\Delta x_{\pic} \Delta y_{\pjc}}
      + \frac{
          \vat{\Delta x \uy}{\pic,\pjp}
        - \vat{\Delta x \uy}{\pic,\pjm}
      }{\Delta x_{\pic} \Delta y_{\pjc}}
   \right) \\
   & = \\
   & - \frac{1}{J_{\pic,\pjc}} \left\{ \diffe{}{\gx} \left( \Delta y \ux \frac{\widehat{T^2}^{\gx}}{2} \right) \right\}_{\pic,\pjc}
     - \frac{1}{J_{\pic,\pjc}} \left\{ \diffe{}{\gy} \left( \Delta x \uy \frac{\widehat{T^2}^{\gy}}{2} \right) \right\}_{\pic,\pjc} \\
   & - \frac{\vat{T^2}{\pic,\pjc}}{2} \left(
      \frac{
          \vat{\ux}{\pip,\pjc}
        - \vat{\ux}{\pim,\pjc}
      }{\Delta x_{\pic}}
      + \frac{
          \vat{\uy}{\pic,\pjp}
        - \vat{\uy}{\pic,\pjm}
      }{\Delta y_{\pjc}}
   \right),

where the last term is null because of the incompressibility constraint.
The first two terms are in conservative forms and thus do not contribute to the change in the total :math:`T^2`.
Note that quadratic quantities are defined as follows:

.. math::

   \vat{\widehat{T^2}^{\gx}}{\pip,\pjc} = \vat{T}{\pic,\pjc} \vat{T}{\pipp,\pjc}, \\
   \vat{\widehat{T^2}^{\gy}}{\pic,\pjp} = \vat{T}{\pic,\pjc} \vat{T}{\pic,\pjpp}.

In summary, I notice that the current advective terms do not contribute to the changes in the total quadratic quantity :math:`H`, i.e. the above scheme is indeed energy-conserving.

***************
Diffusive terms
***************

The local thermal energy dissipation :math:`\dfrac{1}{J} \mst{i}{i} \der{T}{\gx^i} \der{T}{\gx^i}` in the continuous level is given as

.. math::

   \frac{1}{\sqrt{Pr} \sqrt{Ra}} \frac{1}{J} \der{}{\gx^i} \left( \mst{i}{i} T \der{T}{\gx^i} \right)
   -
   \frac{1}{\sqrt{Pr} \sqrt{Ra}} \frac{1}{J} T \der{}{\xi^i} \left( \mst{i}{i} \der{T}{\gx^i} \right).

The first conservative term is discretised as

.. math::

   &
   \frac{1}{\sqrt{Pr} \sqrt{Ra}} \frac{1}{J_{\pic,\pjc}} \left\{
       \vat{\left( \mst{x}{x} T \diffe{T}{\gx} \right)}{\pip,\pjc}
     - \vat{\left( \mst{x}{x} T \diffe{T}{\gx} \right)}{\pim,\pjc}
   \right\} \\
   + &
   \frac{1}{\sqrt{Pr} \sqrt{Ra}} \frac{1}{J_{\pic,\pjc}} \left\{
       \vat{\left( \mst{y}{y} T \diffe{T}{\gy} \right)}{\pic,\pjp}
     - \vat{\left( \mst{y}{y} T \diffe{T}{\gy} \right)}{\pic,\pjm}
   \right\},

while the other term leads

.. math::

   &
   \frac{1}{\sqrt{Pr} \sqrt{Ra}} \frac{1}{J_{\pic,\pjc}} \vat{T}{\pic,\pjc} \left\{
       \vat{\left( \mst{x}{x} \diffe{T}{\gx} \right)}{\pip,\pjc}
     - \vat{\left( \mst{x}{x} \diffe{T}{\gx} \right)}{\pim,\pjc}
   \right\} \\
   + &
   \frac{1}{\sqrt{Pr} \sqrt{Ra}} \frac{1}{J_{\pic,\pjc}} \vat{T}{\pic,\pjc} \left\{
       \vat{\left( \mst{y}{y} \diffe{T}{\gy} \right)}{\pic,\pjp}
     - \vat{\left( \mst{y}{y} \diffe{T}{\gy} \right)}{\pic,\pjm}
   \right\}.

Since

.. math::
   \begin{aligned}
     & + \vat{\left( \dintrpa{T}{\gx} \diffe{T}{\gx} \right)}{\pip,\pjc} - \vat{T}{\pic,\pjc} \vat{\diffe{T}{\gx}}{\pip,\pjc} = \left( + \vat{\dintrpa{T}{\gx}}{\pip,\pjc} - \vat{T}{\pic,\pjc} \right) \vat{\diffe{T}{\gx}}{\pip,\pjc} = \frac{1}{2} \vat{\left( \diffe{T}{\gx} \right)^2}{\pip,\pjc}, \\
     & - \vat{\left( \dintrpa{T}{\gx} \diffe{T}{\gx} \right)}{\pim,\pjc} + \vat{T}{\pic,\pjc} \vat{\diffe{T}{\gx}}{\pim,\pjc} = \left( - \vat{\dintrpa{T}{\gx}}{\pim,\pjc} + \vat{T}{\pic,\pjc} \right) \vat{\diffe{T}{\gx}}{\pim,\pjc} = \frac{1}{2} \vat{\left( \diffe{T}{\gx} \right)^2}{\pim,\pjc}, \\
     & + \vat{\left( \dintrpa{T}{\gy} \diffe{T}{\gy} \right)}{\pic,\pjp} - \vat{T}{\pic,\pjc} \vat{\diffe{T}{\gy}}{\pic,\pjp} = \left( + \vat{\dintrpa{T}{\gy}}{\pic,\pjp} - \vat{T}{\pic,\pjc} \right) \vat{\diffe{T}{\gy}}{\pic,\pjp} = \frac{1}{2} \vat{\left( \diffe{T}{\gy} \right)^2}{\pic,\pjp}, \\
     & - \vat{\left( \dintrpa{T}{\gy} \diffe{T}{\gy} \right)}{\pic,\pjm} + \vat{T}{\pic,\pjc} \vat{\diffe{T}{\gy}}{\pic,\pjm} = \left( - \vat{\dintrpa{T}{\gy}}{\pic,\pjm} + \vat{T}{\pic,\pjc} \right) \vat{\diffe{T}{\gy}}{\pic,\pjm} = \frac{1}{2} \vat{\left( \diffe{T}{\gy} \right)^2}{\pic,\pjm},
   \end{aligned}

the difference of these two terms is

.. math::

   &
   \frac{1}{\sqrt{Pr} \sqrt{Ra}} \frac{1}{J_{\pic,\pjc}} \left[
      \frac{
          \vat{\left\{ \mst{x}{x} \left( \diffe{T}{\gx} \right)^2 \right\}}{\pip,\pjc}
        + \vat{\left\{ \mst{x}{x} \left( \diffe{T}{\gx} \right)^2 \right\}}{\pim,\pjc}
      }{2}
   \right] \\
   +
   &
   \frac{1}{\sqrt{Pr} \sqrt{Ra}} \frac{1}{J_{\pic,\pjc}} \left[
      \frac{
          \vat{\left\{ \mst{y}{y} \left( \diffe{T}{\gy} \right)^2 \right\}}{\pic,\pjp}
        + \vat{\left\{ \mst{y}{y} \left( \diffe{T}{\gy} \right)^2 \right\}}{\pic,\pjm}
      }{2}
   \right] \\
   =
   &
   \frac{1}{\sqrt{Pr} \sqrt{Ra}} \frac{1}{J_{\pic,\pjc}} \left[
      \frac{1}{2} J_{\pip,\pjc} \left( \vat{\frac{\diffe{T}{\gx}}{\Delta x}}{\pip,\pjc} \right)^2
   \right] \\
   + &
   \frac{1}{\sqrt{Pr} \sqrt{Ra}} \frac{1}{J_{\pic,\pjc}} \left[
      \frac{1}{2} J_{\pim,\pjc} \left( \vat{\frac{\diffe{T}{\gx}}{\Delta x}}{\pim,\pjc} \right)^2
   \right] \\
   + &
   \frac{1}{\sqrt{Pr} \sqrt{Ra}} \frac{1}{J_{\pic,\pjc}} \left[
      \frac{1}{2} J_{\pic,\pjp} \left( \vat{\frac{\diffe{T}{\gy}}{\Delta y}}{\pic,\pjp} \right)^2
   \right] \\
   + &
   \frac{1}{\sqrt{Pr} \sqrt{Ra}} \frac{1}{J_{\pic,\pjc}} \left[
      \frac{1}{2} J_{\pic,\pjm} \left( \vat{\frac{\diffe{T}{\gy}}{\Delta y}}{\pic,\pjm} \right)^2
   \right],

which is the local thermal energy dissipation rate.

.. note::

   In the vicinity of the walls, the coefficient :math:`1/2` should be corrected to unity since :math:`T` is defined exactly on the walls.

