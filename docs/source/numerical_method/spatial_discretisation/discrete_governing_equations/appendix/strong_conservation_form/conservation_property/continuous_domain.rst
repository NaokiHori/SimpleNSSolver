#################
Continuous domain
#################

By taking the inner product of the momentum balance and :math:`u_i`, I have

.. math::

   u_i \der{u_i}{t} &
   + u_i \frac{1}{J} \der{}{\xi^j} \left\{ \left( J \der{\xi^j}{x_j} u_j \right) u_i \right\}
   + u_i \frac{1}{J} \der{}{\xi^i} \left( J \der{\xi^i}{x_i} p \right) \\ &
   - u_i \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{1}{J} \der{}{\xi^j} \left( \mst{j}{j} \der{u_i}{\xi^j} \right)
   - u_i f_i
   = 0.

Each term yields the following relation.

*******************
Temporal derivative
*******************

.. math::

   u_i \der{u_i}{t} = \der{}{t} \left( u_i u_i \right) - \der{u_i}{t} u_i

and thus

.. math::

   u_i \der{u_i}{t} = \der{k}{t},

where :math:`k \equiv \frac{1}{2} u_i u_i`.

*********
Advection
*********

.. math::

   \frac{1}{J} \der{}{\xi^i} \left\{ \left( J \der{\xi^i}{x_i} u_i \right) k \right\}.

*****************
Pressure gradient
*****************

.. math::

   \frac{1}{J} \der{}{\xi^i} \left( J \der{\xi^i}{x_i} u_i p \right)
   -
   p \frac{1}{J} \der{}{\xi^i} \left( J \der{\xi^i}{x_i} u_i \right),

where the second term vanishes because of the incompressibility constraint.

*********
Diffusion
*********

.. math::

   - \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{1}{J} \der{}{\xi^j} \left( \mst{j}{j} u_i \der{u_i}{\xi^j} \right)
   + \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{1}{J} \mst{j}{j} \der{u_i}{\xi^j} \der{u_i}{\xi^j}.

**************
Buoyancy force
**************

This term remains as it is.

*******
Summary
*******

In summary, I have the equation of the kinetic energy balance as

.. math::

   \der{k}{t} &
   + \frac{1}{J} \der{}{\xi^i} \left\{ \left( J \der{\xi^i}{x_i} u_i \right) k \right\}
   + \frac{1}{J} \der{}{\xi^i} \left( J \der{\xi^i}{x_i} u_i p \right) \\ &
   - \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{1}{J} \der{}{\xi^j} \left( \mst{j}{j} u_i \der{u_i}{\xi^j} \right)
   + \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{1}{J} \mst{j}{j} \der{u_i}{\xi^j} \der{u_i}{\xi^j}
   - u_i f_i
   = 0.

Similarly, by multiplying the equation of the internal energy by the temperature :math:`T`, I obtain

.. math::

   \der{h}{t} &
   + \frac{1}{J} \der{}{\xi^i} \left\{ \left( J \der{\xi^i}{x_i} u_i \right) h \right\}
   - \frac{1}{\sqrt{Pr} \sqrt{Ra}} \frac{1}{J} \der{}{\xi^j} \left( \mst{j}{j} T \der{T}{\xi^j} \right) \\ &
   + \frac{1}{\sqrt{Pr} \sqrt{Ra}} \frac{1}{J} \mst{j}{j} \der{T}{\xi^j} \der{T}{\xi^j}
   = 0,

where :math:`h \equiv \frac{1}{2} T^2`.

Now I discuss the evolution of the total energies.
The total kinetic energy inside the domain :math:`K` is

.. math::

   \int k dx dy

in the original coordinate system, while

.. math::

   \int k J d\gx d\gy

in the transformed coordinate system.

I obtain the governing equation with respect to :math:`K` by integrating the equation of :math:`k` inside the domain:

.. math::

   \der{K}{t} &
   + \int \der{}{\xi^i} \left\{ \left( J \der{\xi^i}{x_i} u_i \right) k \right\} d\gx d\gy
   + \int \der{}{\xi^i} \left( J \der{\xi^i}{x_i} u_i p \right) d\gx d\gy \\ &
   - \int \frac{\sqrt{Pr}}{\sqrt{Ra}} \der{}{\xi^j} \left( \mst{j}{j} u_i \der{u_i}{\xi^j} \right) d\gx d\gy
   + \int \frac{\sqrt{Pr}}{\sqrt{Ra}} \mst{j}{j} \der{u_i}{\xi^j} \der{u_i}{\xi^j} d\gx d\gy \\ &
   - \int J u_i f_i d\gx d\gy
   = 0.

The second, the third and the fourth terms are conservative and thus the integrated values are determined solely by the boundary fluxes because of the divergence theorem.
For the Rayleigh-Bénard convections considered in this project, periodic boundary conditions or impermeable and no-slip walls (:math:`u_i = 0`) are assumed, and thus these terms vanish.

Thus I have

.. math::

   \der{K}{t}
   + \int \frac{\sqrt{Pr}}{\sqrt{Ra}} \mst{j}{j} \der{u_i}{\xi^j} \der{u_i}{\xi^j} d\gx d\gy
   - \int J u_i f_i d\gx d\gy
   = 0.

When the system is in the statistically-steady state

.. math::

   \ave{\partial K / \partial t}{t} = 0,

I notice that the energy injection by the body force is exactly balanced by the viscous dissipation.

There is an analogue relation for :math:`H = \int h dx dy = \int h J d\gx d\gy`:

.. math::

   \der{H}{t} &
   + \int \der{}{\xi^i} \left\{ \left( J \der{\xi^i}{x_i} u_i \right) h \right\} d\gx d\gy \\ &
   - \int \frac{1}{\sqrt{Pr} \sqrt{Ra}} \der{}{\xi^j} \left( \mst{j}{j} T \der{T}{\xi^j} \right) d\gx d\gy
   + \int \frac{1}{\sqrt{Pr} \sqrt{Ra}} \mst{j}{j} \der{T}{\xi^j} \der{T}{\xi^j} d\gx d\gy
   = 0.

The second term is zero because of the divergence theorem again.
Regarding the third term, it remains here since the thermal energy is continuously injected by the conduction (:math:`\partial T / \partial x`) on the walls.
This injection is exactly dissipated by the last term in a statistical sense.

