
.. _poisson_governing_equation:

##################
Governing equation
##################

For simplicity, in this section, I consider a two-dimensional Poisson equation:

.. math::

   \der{}{x} \der{p}{x}
   +
   \der{}{y} \der{p}{y}
   =
   q,

where

.. math::

   p = p \left( x, y \right),

.. math::

   q = q \left( x, y \right)

and

.. math::

   x \in \left[ 0, l_x \right],

.. math::

   y \in \left[ 0, l_y \right].

The Neumann boundary condition is imposed in the :math:`x` direction:

.. math::

   \vat{\frac{\partial p}{\partial x}}{x =   0}
   =
   \vat{\frac{\partial p}{\partial x}}{x = l_x}
   =
   0,

while the periodicity is assumed in the :math:`y` direction:

.. math::

   \vat{\frac{\partial^n p}{\partial y^n}}{y =   0}
   =
   \vat{\frac{\partial^n p}{\partial y^n}}{y = l_y},

where :math:`n = 0, 1, 2, \cdots`.

In :ref:`the next section <poisson_eigen_decomposition>`, the discretisation of this equation is described and discussed.

