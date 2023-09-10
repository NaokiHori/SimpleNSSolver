##############################
Continuous governing equations
##############################

The non-dimensional governing equations derived :ref:`here <governing_equations>` are, in the strong conservation form, as follows.

.. note::

   :math:`\mst{i}{j} = J \der{\xi^i}{x_k} \der{\xi^j}{x_k}` is the mesh skewness tensor.

****************************
Incompressibility constraint
****************************

.. math::

   \frac{1}{J} \der{}{\xi^i} \left( J \der{\xi^i}{x_i} u_i \right) = 0

********
Momentum
********

.. math::

   \der{u_i}{t}
   + \frac{1}{J} \der{}{\xi^j} \left\{ \left( J \der{\xi^j}{x_j} u_j \right) u_i \right\}
   =
   - \frac{1}{J} \der{}{\xi^i} \left( J \der{\xi^i}{x_i} p \right)
   + \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{1}{J} \der{}{\xi^j} \left( \mst{j}{j} \der{u_i}{\xi^j} \right)
   + T \delta_{ix}

**************
Kinetic energy
**************

.. math::

   \der{k}{t}
   + \frac{1}{J} \der{}{\xi^i} \left\{ \left( J \der{\xi^i}{x_i} u_i \right) \left( k + p \right) \right\}
   =
   \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{1}{J} \der{}{\xi^j} \left( \mst{j}{j} u_i \der{u_i}{\xi^j} \right)
   - \frac{\sqrt{Pr}}{\sqrt{Ra}} \mst{i}{i} \der{u_i}{\xi^j} \der{u_i}{\xi^j}
   + u_i T \delta_{ix}

***************
Internal energy
***************

.. math::

   \der{T}{t}
   + \frac{1}{J} \der{}{\xi^i} \left\{ \left( J \der{\xi^i}{x_i} u_i \right) T \right\}
   = \frac{1}{\sqrt{Pr} \sqrt{Ra}} \frac{1}{J} \der{}{\xi^i} \left( \mst{i}{i} \der{T}{\xi^i} \right)

**************************
Thermal quadratic quantity
**************************

.. math::

   \der{h}{t}
   + \frac{1}{J} \der{}{\xi^i} \left\{ \left( J \der{\xi^i}{x_i} u_i \right) h \right\}
   =
   \frac{1}{\sqrt{Pr} \sqrt{Ra}} \frac{1}{J} \der{}{\xi^i} \left( \mst{i}{i} T \der{T}{\xi^i} \right)
   - \frac{1}{\sqrt{Pr} \sqrt{Ra}} \mst{i}{i} \der{T}{\xi^j} \der{T}{\xi^j}

