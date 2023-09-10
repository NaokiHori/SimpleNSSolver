
.. include:: /references.txt

.. _implicit_treatment:

##################
Implicit treatment
##################

****************************
Time scales and restrictions
****************************

Numerically, I should not allow any information to propagate longer distance than the grid size in one time step.
To use a larger time step size :math:`\Delta t`, the information should be treated implicitly which is capable of capturing infinite propagation speed.
There are following three time scales in the governing equations which are considered in this project, which give the corresponding three restrictions on the time step size :math:`\Delta t`.

#. Advective terms

   One comes from the advective terms in the governing equations, where the information should not travel longer distance than the local grid size.
   This constraint yields

   .. math::

      \Delta t_{adv}
      <
      \Delta x / \ux,

   where :math:`\Delta x` and :math:`\ux` are the local grid size and the velocity, respectively.
   Note that, although written only in :math:`x` direction here, the same condition is applied to the other directions.

   I see that the time step size should be reduced as the resolution gets finer, and the refinement speed is proportional to the grid size:

   .. math::

      \Delta t_{adv}
      =
      C \Delta x,

   where :math:`C` is called Courant number.

   Since the advective terms are non-linear, it is impractical to treat them implicitly, and I cannot elimitate this constraint.

#. Diffusive terms

   The other one is imposed by the diffusive terms in the governing equations, where the time scale yields

   .. math::

      & \Delta t_{dif, fluid}       \propto \frac{\sqrt{Ra}}{\sqrt{Pr}} \Delta x^2, \\
      & \Delta t_{dif, temperature} \propto \sqrt{Pr} \sqrt{Ra}         \Delta x^2,

   in non-dimensional form.
   Note that the kinematic viscosity has the unit of :math:`\left[ L^2 T^{-1} \right]`, giving the characteristic time scale above.

   Now I see that the time step size should be reduced quadratically:

   .. math::

      \Delta t_{dif}
      =
      F \Delta x^2,

   where :math:`F` is called Fourier number.

   This criterion is often very severe and easily makes computations impractical, especially for wall-bounded turbulent flows where the wall-normal grid sizes should be extremely small close to the walls to resolve the boundary layers.

   Since the diffusive terms are linear, this restriction can be eliminated by treating it implicitly, which is the central focus of this section.

#. Pressure-gradient terms

   Incompressible liquids is the limited condition whose speed of sound is infinity.
   To avoid :math:`\Delta t = 0`, the pressure field should always be treated implicitly.

   .. seealso::

      In other words, for incompressible flows, the pressure field is determined to satisfy the incompressibility, which is the central idea of :ref:`the SMAC method <smac_method>`.

*************************************************
Approximate factorisation for the diffusive terms
*************************************************

To treat the diffusive terms implicitly, I need to solve the following N-dimensional Helmholtz equations:

.. math::

   \left( 1 - \frac{\Delta t}{2} \frac{\sqrt{Pr}}{\sqrt{Ra}} \dder{}{x_j} \dder{}{x_j} \right) u_i^{n+1}
   =
   -
   \dder{p^n}{x_i} \Delta t
   -
   u_j^n \dder{u_i^n}{x_j} \Delta t
   +
   \left( 1 + \frac{\Delta t}{2} \frac{\sqrt{Pr}}{\sqrt{Ra}} \dder{}{x_j} \dder{}{x_j} \right) u_i^n

for the momentum field, and

.. math::

   \left( 1 - \frac{\Delta t}{2} \frac{1}{\sqrt{Pr} \sqrt{Ra}} \dder{}{x_j} \dder{}{x_j} \right) T^{n+1}
   =
   -
   u_j^n \dder{T^n}{x_j} \Delta t
   +
   \left( 1 + \frac{\Delta t}{2} \frac{1}{\sqrt{Pr} \sqrt{Ra}} \dder{}{x_j} \dder{}{x_j} \right) T^n

for the temperature field.

In this project, I can simplify these equations as follows.

First, I re-write the equations as

.. math::

   \left( 1 - \frac{\Delta t}{2} \frac{\sqrt{Pr}}{\sqrt{Ra}} \dder{}{x_j} \dder{}{x_j} \right) \Delta u_i
   =
   -
   \dder{p^n}{x_i} \Delta t
   -
   u_j^n \dder{u_i^n}{x_j} \Delta t
   +
   \Delta t \frac{\sqrt{Pr}}{\sqrt{Ra}} \dder{}{x_j} \dder{}{x_j} u_i^n,

.. math::

   \left( 1 - \frac{\Delta t}{2} \frac{1}{\sqrt{Pr} \sqrt{Ra}} \dder{}{x_j} \dder{}{x_j} \right) \Delta T
   =
   -
   u_j^n \dder{T^n}{x_j} \Delta t
   +
   \Delta t \frac{1}{\sqrt{Pr} \sqrt{Ra}} \dder{}{x_j} \dder{}{x_j} T^n,

where

.. math::

   \Delta u_i
   \equiv
   u_i^{n+1}
   -
   u_i^{n  },

.. math::

   \Delta T
   \equiv
   T^{n+1}
   -
   T^{n  }.

Then I approximate the left-hand-side terms as

.. math::

   \left( 1 - \frac{\Delta t}{2} \frac{\sqrt{Pr}}{\sqrt{Ra}} \dder{}{x_j} \dder{}{x_j} \right) \Delta u_i
   &
   =
   \left(
      1
      -
      \frac{\Delta t}{2} \frac{\sqrt{Pr}}{\sqrt{Ra}} \dder{}{x} \dder{}{x}
      -
      \frac{\Delta t}{2} \frac{\sqrt{Pr}}{\sqrt{Ra}} \dder{}{y} \dder{}{y}
   \right) \Delta u_i \\
   &
   \approx
   \left( 1 - \frac{\Delta t}{2} \frac{\sqrt{Pr}}{\sqrt{Ra}} \dder{}{x} \dder{}{x} \right)
   \left( 1 - \frac{\Delta t}{2} \frac{\sqrt{Pr}}{\sqrt{Ra}} \dder{}{y} \dder{}{y} \right)
   \Delta u_i,

.. math::

   \left( 1 - \frac{\Delta t}{2} \frac{1}{\sqrt{Pr} \sqrt{Ra}} \dder{}{x_j} \dder{}{x_j} \right) \Delta T
   &
   =
   \left(
      1
      -
      \frac{\Delta t}{2} \frac{1}{\sqrt{Pr} \sqrt{Ra}} \dder{}{x} \dder{}{x}
      -
      \frac{\Delta t}{2} \frac{1}{\sqrt{Pr} \sqrt{Ra}} \dder{}{y} \dder{}{y}
   \right) \Delta T \\
   &
   \approx
   \left( 1 - \frac{\Delta t}{2} \frac{1}{\sqrt{Pr} \sqrt{Ra}} \dder{}{x} \dder{}{x} \right)
   \left( 1 - \frac{\Delta t}{2} \frac{1}{\sqrt{Pr} \sqrt{Ra}} \dder{}{y} \dder{}{y} \right)
   \Delta T.

By assuming

.. math::

   \Delta u_i
   \sim
   \Delta t,
   \Delta T
   \sim
   \Delta t,

I notice that the splitting errors

.. math::

   \frac{\Delta t}{2} \frac{\sqrt{Pr}}{\sqrt{Ra}} \dder{}{x} \dder{}{x}
   \times
   \frac{\Delta t}{2} \frac{\sqrt{Pr}}{\sqrt{Ra}} \dder{}{y} \dder{}{y}
   \times
   \Delta u_i,

.. math::

   \frac{\Delta t}{2} \frac{1}{\sqrt{Pr} \sqrt{Ra}} \dder{}{x} \dder{}{x}
   \times
   \frac{\Delta t}{2} \frac{1}{\sqrt{Pr} \sqrt{Ra}} \dder{}{y} \dder{}{y}
   \times
   \Delta T,

are :math:`\sim \Delta t^3`, which is comparable to the dominant error of :ref:`the used Runge-Kutta scheme <time_marchers>`.

In summary, the N-dimensional Helmholtz equations are approximated by the linear system for each direction and much easier to solve.
This splitting method is called approximate factorisation (|DUKOWICZ1992|).
Obviously this technique is not accepted if

   * I solve the equation with respect to the variable itself (i.e. :math:`u_i` instead of :math:`\Delta u_i`)

   * I use higher-order scheme to integrate the equation in time.

Here, I need to solve :ref:`a linear system <linear_system>` to find :math:`\left( \cdots \right)^{-1}`.
Since they are independent linear systems in each direction and I adopt the second-order-accurate central-difference scheme in space, I can solve them by :ref:`the tri-diagonal matrix algorithm <tdm>`.

.. note::

   * Overhead

      In this project, the implicit treatment in the :math:`x` (wall-normal) direction can be easily achieved, whose cost is up to a few percent.
      The implicit treatments in the other directions, however, require certain amount of MPI communication, whose overhead can be more than :math:`100` percent.

   * Monotonicity

      Although the Crank-Nicolson scheme can eliminate the stability restriction, monotonicity (i.e. temperature is bounded between the two boundary values) is not guaranteed.
      Unfortunately, in order to guarantee the monotonicity, :math:`\Delta t \propto \Delta x^2` should be satisfied again (see e.g. |HORVATH2000|).
      Although this restriction disappears by adopting the Euler backward scheme (at the expense of the temporal accuracy), this fact is neglected in this project for now.
      To minimises this issue, using a sufficiently fine spatial resolution to resolve the highest frequency is important.

