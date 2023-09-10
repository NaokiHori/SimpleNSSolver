#######################
Discrete kinetic energy
#######################

I consider to multiply the momentum balance in the :math:`x` direction by :math:`\vat{\ux}{\xic,\xjc}`, yielding

.. math::

   \vat{\ux}{\xic,\xjc} \der{\vat{\ux}{\xic,\xjc}}{t}
   =
   & - \vat{\ux}{\xic,\xjc} \frac{
       \vat{\dintrpa{\Delta y \ux}{\gx}}{\xip,\xjc} \vat{\dintrpa{\ux}{\gx}}{\xip,\xjc}
     - \vat{\dintrpa{\Delta y \ux}{\gx}}{\xim,\xjc} \vat{\dintrpa{\ux}{\gx}}{\xim,\xjc}
   }{\Delta x_{\xic} \Delta y_{\xjc}} \\
   & - \vat{\ux}{\xic,\xjc} \frac{
       \vat{\dintrpa{\Delta x \uy}{\gx}}{\xic,\xjp} \vat{\dintrpa{\ux}{\gy}}{\xic,\xjp}
     - \vat{\dintrpa{\Delta x \uy}{\gx}}{\xic,\xjm} \vat{\dintrpa{\ux}{\gy}}{\xic,\xjm}
   }{\Delta x_{\xic} \Delta y_{\xjc}} \\
   & - \vat{\ux}{\xic,\xjc} \frac{
       \vat{\Delta y p}{\xip,\xjc}
     - \vat{\Delta y p}{\xim,\xjc}
   }{\Delta x_{\xic} \Delta y_{\xjc}} \\
   & + \vat{\ux}{\xic,\xjc} \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{
       \vat{\frac{\Delta y}{\Delta x} \diffe{\ux}{\gx}}{\xip,\xjc}
     - \vat{\frac{\Delta y}{\Delta x} \diffe{\ux}{\gx}}{\xim,\xjc}
   }{\Delta x_{\xic} \Delta y_{\xjc}} \\
   & + \vat{u}{\xic,\xjc} \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{
       \vat{\frac{\Delta x}{\Delta y} \diffe{\ux}{\gy}}{\xic,\xjp}
     - \vat{\frac{\Delta x}{\Delta y} \diffe{\ux}{\gy}}{\xic,\xjm}
   }{\Delta x_{\xic} \Delta y_{\xjc}} \\
   & + \vat{\ux}{\xic,\xjc} \vat{\dintrpa{T}{\gx}}{\xic,\xjc}.

Similarly, the momentum balance in the :math:`y` direction is multiplied by :math:`\vat{\uy}{\yic,\yjc}`:

.. math::

   \vat{\uy}{\yic,\yjc} \der{\vat{\uy}{\yic,\yjc}}{t}
   =
   & - \vat{\uy}{\yic,\yjc} \frac{
       \vat{\dintrpa{\Delta y \ux}{\gy}}{\yip,\yjc} \vat{\dintrpa{\uy}{\gx}}{\yip,\yjc}
     - \vat{\dintrpa{\Delta y \ux}{\gy}}{\yim,\yjc} \vat{\dintrpa{\uy}{\gx}}{\yim,\yjc}
   }{\Delta x_{\yic} \Delta y_{\yjc}} \\
   & - \vat{\uy}{\yic,\yjc} \frac{
       \vat{\dintrpa{\Delta x \uy}{\gy}}{\yic,\yjp} \vat{\dintrpa{\uy}{\gy}}{\yic,\yjp}
     - \vat{\dintrpa{\Delta x \uy}{\gy}}{\yic,\yjm} \vat{\dintrpa{\uy}{\gy}}{\yic,\yjm}
   }{\Delta x_{\yic} \Delta y_{\yjc}} \\
   & - \vat{\uy}{\yic,\yjc} \frac{
       \vat{\Delta x p}{\yic,\yjp}
     - \vat{\Delta x p}{\yic,\yjm}
   }{\Delta x_{\yic} \Delta y_{\yjc}} \\
   & + \vat{\uy}{\yic,\yjc} \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{
       \vat{\frac{\Delta y}{\Delta x} \diffe{\uy}{\gx}}{\yip,\yjc}
     - \vat{\frac{\Delta y}{\Delta x} \diffe{\uy}{\gx}}{\yim,\yjc}
   }{\Delta x_{\yic} \Delta y_{\yjc}} \\
   & + \vat{\uy}{\yic,\yjc} \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{
       \vat{\frac{\Delta x}{\Delta y} \diffe{\uy}{\gy}}{\yic,\yjp}
     - \vat{\frac{\Delta x}{\Delta y} \diffe{\uy}{\gy}}{\yic,\yjm}
   }{\Delta x_{\yic} \Delta y_{\yjc}}.

The left-hand-side terms give the change of :math:`1/2 \ux^2` and :math:`1/2 \uy^2` in time, which are clearly relevant to the discrete kinetic energy.

***************
Advective terms
***************

==============
Local property
==============

Since the velocity components are located at different positions, the discrete kinetic energy cannot be defined uniquely.
Here, I consider the squared velocities (:math:`\ux^2` and :math:`\uy^2`) separately at positions where they are defined.

In the :math:`x` direction, the advection terms at :math:`\left( \xic, \xjc \right)` lead to

.. math::

   & - \frac{
       \vat{\dintrpa{\Delta y \ux}{\gx}}{\xip,\xjc} \frac{\vat{\ux}{\xipp,\xjc } \vat{\ux}{\xic ,\xjc }}{2}
     - \vat{\dintrpa{\Delta y \ux}{\gx}}{\xim,\xjc} \frac{\vat{\ux}{\xic ,\xjc } \vat{\ux}{\ximm,\xjc }}{2}
   }{\Delta x_{\xic} \Delta y_{\xjc}} \\
   & - \frac{
       \vat{\dintrpa{\Delta x \uy}{\gx}}{\xic,\xjp} \frac{\vat{\ux}{\xic ,\xjpp} \vat{\ux}{\xic ,\xjc }}{2}
     - \vat{\dintrpa{\Delta x \uy}{\gx}}{\xic,\xjm} \frac{\vat{\ux}{\xic ,\xjc } \vat{\ux}{\xic ,\xjmm}}{2}
   }{\Delta x_{\xic} \Delta y_{\xjc}} \\
   & - \frac{\vat{\ux^2}{\xic,\xjc}}{2} \frac{
        \vat{\dintrpa{\Delta y \ux}{\gx}}{\xip,\xjc}
      - \vat{\dintrpa{\Delta y \ux}{\gx}}{\xim,\xjc}
      + \vat{\dintrpa{\Delta x \uy}{\gx}}{\xic,\xjp}
      - \vat{\dintrpa{\Delta x \uy}{\gx}}{\xic,\xjm}
   }{\Delta x_{\xic} \Delta y_{\xjc}} \\
   & = \\
   & -\frac{1}{J_{\xic,\xjc}} \left\{ \diffe{}{\gx} \left(
      \dintrpa{\Delta y \ux}{\gx} \frac{{\widehat{\ux^2}}^{\gx}}{2}
   \right) \right\}_{\xic,\xjc}
   - \frac{1}{J_{\xic,\xjc}} \left\{ \diffe{}{\gy} \left(
      \dintrpa{\Delta x \uy}{\gx} \frac{{\widehat{\ux^2}}^{\gy}}{2}
   \right) \right\}_{\xic,\xjc} \\
   & - \frac{\vat{u^2}{\xic,\xjc}}{2} \frac{1}{J_{\xic,\xjc}} \left\{
      \frac{J_{\xim,\xjc}}{2} \left(
        \frac{
            \vat{\ux}{\xic, \xjc}
          - \vat{\ux}{\ximm,\xjc}
        }{\Delta x_{\xim}}
        + \frac{
            \vat{\uy}{\xim,\xjp}
          - \vat{\uy}{\xim,\xjm}
        }{\Delta y_{\xjc}}
      \right)
   \right\} \\
   & - \frac{\vat{u^2}{\xic,\xjc}}{2} \frac{1}{J_{\xic,\xjc}} \left\{
      \frac{J_{\xip,\xjc}}{2} \left(
        \frac{
            \vat{\ux}{\xipp,\xjc}
          - \vat{\ux}{\xic, \xjc}
        }{\Delta x_{\xip}}
        + \frac{
            \vat{\uy}{\xip,\xjp}
          - \vat{\uy}{\xip,\xjm}
        }{\Delta y_{\xjc}}
      \right)
   \right\},

while in the :math:`y` direction at :math:`\left( \yic, \yjc \right)`, I have

.. math::

   & - \frac{
       \vat{\dintrpa{\Delta y \ux}{\gy}}{\yip,\yjc} \frac{\vat{\uy}{\yipp,\yjc } \vat{\uy}{\yic, \yjc }}{2}
     - \vat{\dintrpa{\Delta y \ux}{\gy}}{\yim,\yjc} \frac{\vat{\uy}{\yic ,\yjc } \vat{\uy}{\yimm,\yjc }}{2}
   }{\Delta x_{\yic} \Delta y_{\yjc}} \\
   & - \frac{
       \vat{\dintrpa{\Delta x \uy}{\gy}}{\yic,\yjp} \frac{\vat{\uy}{\yic ,\yjpp} \vat{\uy}{\yic ,\yjc }}{2}
     - \vat{\dintrpa{\Delta x \uy}{\gy}}{\yic,\yjm} \frac{\vat{\uy}{\yic ,\yjc } \vat{\uy}{\yic ,\yjmm}}{2}
   }{\Delta x_{\yic} \Delta y_{\yjc}} \\
   & - \frac{\vat{v^2}{\yic,\yjc}}{2} \frac{
        \vat{\dintrpa{\Delta y \ux}{\gy}}{\yip,\yjc}
      - \vat{\dintrpa{\Delta y \ux}{\gy}}{\yim,\yjc}
      + \vat{\dintrpa{\Delta x \uy}{\gy}}{\yic,\yjp}
      - \vat{\dintrpa{\Delta x \uy}{\gy}}{\yic,\yjm}
   }{\Delta x_{\yic} \Delta y_{\yjc}} \\
   & = \\
   & -\frac{1}{J_{\yic,\yjc}} \left\{ \diffe{}{\gx} \left(
      \dintrpa{\Delta y \ux}{\gx} \frac{{\widehat{\uy^2}}^{\gx}}{2}
   \right) \right\}_{\yic,\yjc}
   - \frac{1}{J_{\yic,\yjc}} \left\{ \diffe{}{\gy} \left(
      \dintrpa{\Delta x \uy}{\gx} \frac{{\widehat{\uy^2}}^{\gy}}{2}
   \right) \right\}_{\yic,\yjc} \\
   & - \frac{\vat{\uy^2}{\yic,\yjc}}{2} \frac{1}{J_{\yic,\yjc}} \left\{
      \frac{J_{\yic,\yjm}}{2} \left(
        \frac{
            \vat{\ux}{\yim,\yjm}
          - \vat{\ux}{\yip,\yjm}
        }{\Delta x_{\yic}}
        + \frac{
            \vat{\uy}{\yic,\yjmm}
          - \vat{\uy}{\yic,\yjc }
        }{\Delta y_{\yjm}}
      \right)
   \right\} \\
   & - \frac{\vat{\uy^2}{\yic,\yjc}}{2} \frac{1}{J_{\yic,\yjc}} \left\{
      \frac{J_{\yic,\yjp}}{2} \left(
        \frac{
            \vat{\ux}{\yim,\yjp}
          - \vat{\ux}{\yip,\yjp}
        }{\Delta x_{\yic}}
        + \frac{
            \vat{\uy}{\yic,\yjpp}
          - \vat{\uy}{\yic,\yjc }
        }{\Delta y_{\yjp}}
      \right)
   \right\}.

Note that the last terms are the volume average of the incompressibility constraints and thus zero.

Here new symbols are introduced (quadratic quantities, which are the pseudo kinetic energies):

.. math::

   & \vat{\widehat{\ux^2}^{\gx}}{\xip,\xjc} \equiv \vat{\ux}{\xic,\xjc} \vat{\ux}{\xipp,\xjc }, \\
   & \vat{\widehat{\ux^2}^{\gy}}{\xic,\xjp} \equiv \vat{\ux}{\xic,\xjc} \vat{\ux}{\xic ,\xjpp},

in the :math:`x` direction, while

.. math::

   & \vat{\widehat{\uy^2}^{\gx}}{\yip,\yjc} \equiv \vat{\uy}{\yic,\yjc} \vat{v}{\yipp,\yjc }, \\
   & \vat{\widehat{\uy^2}^{\gy}}{\yic,\yjp} \equiv \vat{\uy}{\yic,\yjc} \vat{v}{\yic ,\yjpp},

in the :math:`y` direction.

===============
Global property
===============

Now I consider to integrate the above two advective terms in the whole domain.

.. note::

   The integral of a quantity :math:`q`

   .. math::

      \int q dx dy = \int q J d\gx d\gy

   is discretised as

   .. math::

      \sum q_{i,j} \Delta x_i \Delta y_j = \sum q_{i,j} J_{i,j} \Delta \gx_i \Delta \gy_j = \sum q_{i,j} J_{i,j} \,\,\, \left( \because \Delta \gx_i \equiv \Delta \gy_j \equiv 1 \right)

   in the general coordinate system.

:math:`\der{K}{t}`, which is the evolution of the total (global) discrete kinetic energy :math:`K`, is given by the sum of

.. math::

   \sum_{\forall \ux \text{positions}, \left( \xic, \xjc \right)} \der{}{t} \left( \frac{1}{2} J \ux^2 \right)_{\xic,\xjc}
   \equiv
   \sum_{\forall \ux \text{positions}, \left( \xic, \xjc \right)}
   \left[
      \begin{aligned}
         & - \left\{ \diffe{}{\gx} \left(
            \dintrpa{\Delta y \ux}{\gx} \frac{{\widehat{\ux^2}}^{\gx}}{2}
         \right) \right\}_{\xic,\xjc} \\
         & - \left\{ \diffe{}{\gy} \left(
            \dintrpa{\Delta x \uy}{\gx} \frac{{\widehat{\ux^2}}^{\gy}}{2}
         \right) \right\}_{\xic,\xjc} \\
         & - \vat{\ux}{\xic,\xjc} \Delta y_{\xjc} \left(
             \vat{p}{\xip,\xjc}
           - \vat{p}{\xim,\xjc}
         \right) \\
      \end{aligned}
   \right]

and

.. math::

   \sum_{\forall \uy \text{positions}, \left( \yic, \yjc \right)} \der{}{t} \left( \frac{1}{2} J \uy^2 \right)_{\yic,\yjc}
   \equiv
   \sum_{\forall \uy \text{positions}, \left( \yic, \yjc \right)}
   \left[
      \begin{aligned}
         & - \left\{ \diffe{}{\gx} \left(
            \dintrpa{\Delta y \ux}{\gx} \frac{{\widehat{\uy^2}}^{\gx}}{2}
         \right) \right\}_{\yic,\yjc} \\
         & - \left\{ \diffe{}{\gy} \left(
            \dintrpa{\Delta x \uy}{\gx} \frac{{\widehat{\uy^2}}^{\gy}}{2}
         \right) \right\}_{\yic,\yjc} \\
         & - \vat{\uy}{\yic,\yjc} \Delta x_{\yic} \left(
             \vat{p}{\yic,\yjp}
           - \vat{p}{\yic,\yjm}
         \right) \\
      \end{aligned}
   \right].

In each summation, the first two terms vanish since they are conservative.
Regarding the pressure gradient terms, they can be re-ordered as

.. math::

   \sum_{\forall p \text{positions}, \left( \pic, \pjc \right)} \vat{p}{\pic,\pjc} \Delta y_j \left( \vat{\ux}{\pip,\pjc} - \vat{\ux}{\pim,\pjc} \right)
   +
   \sum_{\forall p \text{positions}, \left( \pic, \pjc \right)} \vat{p}{\pic,\pjc} \Delta x_i \left( \vat{\uy}{\pic,\pjp} - \vat{\uy}{\pic,\pjm} \right),

which is

.. math::

   \sum_{\forall p \text{positions}, \left( \pic, \pjc \right)} \vat{p J}{\pic,\pjc} \left(
       \frac{\vat{\ux}{\pip,\pjc} - \vat{\ux}{\pim,\pjc}}{\Delta x_i}
     + \frac{\vat{\uy}{\pic,\pjp} - \vat{\uy}{\pic,\pjm}}{\Delta y_j}
   \right),

and thus null because of the incompressibility constraint again.

In summary, I notice that the current advective and pressure gradient terms do not contribute to the changes in the total quadratic quantity :math:`K`, i.e. the above scheme is indeed energy-conserving.

***************
Diffusive terms
***************

The diffusive terms play two roles:

   * transporting the kinetic energy,

   * dissipating the kinetic energy.

The latter is important when computing the kinetic energy dissipation (and the Nusselt number).
To obtain a consistent discretisation, I try to mimic the relation in the continuous domain: the viscous dissipation arising from the momentum equation in the :math:`x` direction is given as

.. math::

   \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{1}{J} \mst{j}{j} \der{\ux}{\xi^j} \der{\ux}{\xi^j}
   = \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{1}{J} \der{}{\xi^j} \left( \mst{j}{j} \ux \der{\ux}{\xi^j} \right)
   - \ux \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{1}{J} \der{}{\xi^j} \left( \mst{j}{j} \der{\ux}{\xi^j} \right),

which is directly discretised:

.. math::

   & \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{1}{J_{\xic,\xjc}} \left\{
       \vat{\left( \mst{x}{x} \dintrpa{\ux}{\gx} \right) \diffe{\ux}{\gx}}{\xip,\xjc}
     - \vat{\left( \mst{x}{x} \dintrpa{\ux}{\gx} \right) \diffe{\ux}{\gx}}{\xim,\xjc}
   \right\} \\
   +
   & \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{1}{J_{\xic,\xjc}} \left\{
       \vat{\left( \mst{y}{y} \dintrpa{\ux}{\gy} \right) \diffe{\ux}{\gy}}{\xic,\xjp}
     - \vat{\left( \mst{y}{y} \dintrpa{\ux}{\gy} \right) \diffe{\ux}{\gy}}{\xic,\xjm}
   \right\} \\
   -
   & \vat{\ux}{\xic,\xjc} \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{1}{J_{\xic,\xjc}} \left\{
       \vat{\left( \mst{x}{x} \diffe{\ux}{\gx} \right)}{\xip,\xjc}
     - \vat{\left( \mst{x}{x} \diffe{\ux}{\gx} \right)}{\xim,\xjc}
   \right\} \\
   -
   & \vat{\ux}{\xic,\xjc} \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{1}{J_{\xic,\xjc}} \left\{
       \vat{\left( \mst{y}{y} \diffe{\ux}{\gy} \right)}{\xic,\xjp}
     - \vat{\left( \mst{y}{y} \diffe{\ux}{\gy} \right)}{\xic,\xjm}
   \right\},

which is

.. math::

   & + \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{1}{J_{\xic,\xjc}} \left( \vat{\dintrpa{\ux}{\gx}}{\xip,\xjc} - \vat{\ux}{\xic,\xjc} \right) \vat{\left( \mst{x}{x} \diffe{\ux}{\gx} \right)}{\xip,\xjc} \\
   & - \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{1}{J_{\xic,\xjc}} \left( \vat{\dintrpa{\ux}{\gx}}{\xim,\xjc} - \vat{\ux}{\xic,\xjc} \right) \vat{\left( \mst{x}{x} \diffe{\ux}{\gx} \right)}{\xim,\xjc} \\
   & + \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{1}{J_{\xic,\xjc}} \left( \vat{\dintrpa{\ux}{\gy}}{\xic,\xjp} - \vat{\ux}{\xic,\xjc} \right) \vat{\left( \mst{y}{y} \diffe{\ux}{\gy} \right)}{\xic,\xjp} \\
   & - \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{1}{J_{\xic,\xjc}} \left( \vat{\dintrpa{\ux}{\gy}}{\xic,\xjm} - \vat{\ux}{\xic,\xjc} \right) \vat{\left( \mst{y}{y} \diffe{\ux}{\gy} \right)}{\xic,\xjm}.

Since I have

.. math::
   \begin{aligned}
     & \vat{\dintrpa{\ux}{\gx}}{\xip,\xjc} - \vat{\ux}{\xic,\xjc}
       = \frac{\vat{\ux}{\xipp,\xjc} + \vat{\ux}{\xic ,\xjc}}{2} - \vat{\ux}{\xic,\xjc}
       = + \frac{\vat{\ux}{\xipp,\xjc} - \vat{\ux}{\xic ,\xjc}}{2}
       = + \frac{1}{2} \vat{\diffe{\ux}{\gx}}{\xip,\xjc} \\
     & \vat{\dintrpa{\ux}{\gx}}{\xim,\xjc} - \vat{\ux}{\xic,\xjc}
       = \frac{\vat{\ux}{\xic ,\xjc} + \vat{\ux}{\ximm,\xjc}}{2} - \vat{\ux}{\xic,\xjc}
       = - \frac{\vat{\ux}{\xic ,\xjc} - \vat{\ux}{\ximm,\xjc}}{2}
       = - \frac{1}{2} \vat{\diffe{\ux}{\gx}}{\xim,\xjc} \\
     & \vat{\dintrpa{\ux}{\gy}}{\xic,\xjp} - \vat{\ux}{\xic,\xjc}
       = \frac{\vat{\ux}{\xic,\xjpp} + \vat{\ux}{\xic,\xjc }}{2} - \vat{\ux}{\xic,\xjc}
       = + \frac{\vat{\ux}{\xic,\xjpp} - \vat{\ux}{\xic,\xjc }}{2}
       = + \frac{1}{2} \vat{\diffe{\ux}{\gy}}{\xic,\xjp} \\
     & \vat{\dintrpa{\ux}{\gy}}{\xic,\xjm} - \vat{\ux}{\xic,\xjc}
       = \frac{\vat{\ux}{\xic,\xjc } + \vat{\ux}{\xic,\xjmm}}{2} - \vat{\ux}{\xic,\xjc}
       = - \frac{\vat{\ux}{\xic,\xjc } - \vat{\ux}{\xic,\xjmm}}{2}
       = - \frac{1}{2} \vat{\diffe{\ux}{\gy}}{\xic,\xjm}
   \end{aligned}

in the bulk, I finally notice that the viscous dissipation arising from the momentum equation in the :math:`x` direction leads to

.. math::

   & \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{1}{J_{\xic,\xjc}} \left[
      \frac{
          \vat{\left\{ \mst{x}{x} \left( \diffe{\ux}{\gx} \right)^2 \right\}}{\xip,\xjc}
        + \vat{\left\{ \mst{x}{x} \left( \diffe{\ux}{\gx} \right)^2 \right\}}{\xim,\xjc}
      }{2}
      + \frac{
          \vat{\left\{ \mst{y}{y} \left( \diffe{\ux}{\gy} \right)^2 \right\}}{\xic,\xjp}
        + \vat{\left\{ \mst{y}{y} \left( \diffe{\ux}{\gy} \right)^2 \right\}}{\xic,\xjm}
      }{2}
   \right] \\
   =
   & \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{1}{J_{\xic,\xjc}} \left\{
      \frac{1}{2} J_{\xip,\xjc} \left( \vat{\frac{\diffe{\ux}{\gx}}{\Delta x}}{\xip,\xjc} \right)^2
   \right\} \\
   +
   & \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{1}{J_{\xic,\xjc}} \left\{
      \frac{1}{2} J_{\xim,\xjc} \left( \vat{\frac{\diffe{\ux}{\gx}}{\Delta x}}{\xim,\xjc} \right)^2
   \right\} \\
   +
   & \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{1}{J_{\xic,\xjc}} \left\{
      \frac{1}{2} J_{\xic,\xjp} \left( \vat{\frac{\diffe{\ux}{\gy}}{\Delta y}}{\xic,\xjp} \right)^2
   \right\} \\
   +
   & \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{1}{J_{\xic,\xjc}} \left\{
      \frac{1}{2} J_{\xic,\xjm} \left( \vat{\frac{\diffe{\ux}{\gy}}{\Delta y}}{\xic,\xjm} \right)^2
   \right\}

at :math:`\left( \xic, \xjc \right).`

Following a similar manner, by using the relation

.. math::

   \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{1}{J} \mst{j}{j} \der{\uy}{\xi^j} \der{\uy}{\xi^j}
   = \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{1}{J} \der{}{\xi^j} \left( \mst{j}{j} \uy \der{\uy}{\xi^j} \right)
   - \uy \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{1}{J} \der{}{\xi^j} \left( \mst{j}{j} \der{\uy}{\xi^j} \right),

I can derive the viscous dissipation coming from the momentum equation in the :math:`y` direction as

.. math::

   & \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{1}{J_{\yic,\yjc}} \left[
      \frac{
          \vat{\left\{ \mst{x}{x} \left( \diffe{\uy}{\gx} \right)^2 \right\}}{\yip,\yjc}
        + \vat{\left\{ \mst{x}{x} \left( \diffe{\uy}{\gx} \right)^2 \right\}}{\yim,\yjc}
      }{2}
      + \frac{
          \vat{\left\{ \mst{y}{y} \left( \diffe{\uy}{\gy} \right)^2 \right\}}{\yic,\yjp}
        + \vat{\left\{ \mst{y}{y} \left( \diffe{\uy}{\gy} \right)^2 \right\}}{\yic,\yjm}
      }{2}
   \right] \\
   =
   & \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{1}{J_{\yic,\yjc}} \left\{
      \frac{1}{2} J_{\yip,\yjc} \left( \vat{\frac{\diffe{\uy}{\gx}}{\Delta x}}{\yip,\yjc} \right)^2
   \right\} \\
   +
   & \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{1}{J_{\yic,\yjc}} \left\{
      \frac{1}{2} J_{\yim,\yjc} \left( \vat{\frac{\diffe{\uy}{\gx}}{\Delta x}}{\yim,\yjc} \right)^2
   \right\} \\
   +
   & \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{1}{J_{\yic,\yjc}} \left\{
      \frac{1}{2} J_{\yic,\yjp} \left( \vat{\frac{\diffe{\uy}{\gy}}{\Delta y}}{\yic,\yjp} \right)^2
   \right\} \\
   +
   & \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{1}{J_{\yic,\yjc}} \left\{
      \frac{1}{2} J_{\yic,\yjm} \left( \vat{\frac{\diffe{\uy}{\gy}}{\Delta y}}{\yic,\yjm} \right)^2
   \right\}

at :math:`\left( \yic, \yjc \right).`

.. note::

   For the :math:`y` velocities, in the vicinity of the walls, the coefficient :math:`1/2` should be corrected to :math:`1` since :math:`\uy` is defined exactly on the walls.

