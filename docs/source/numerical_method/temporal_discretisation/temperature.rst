
.. _temperature_integration:

####################################
Integrating internal energy equation
####################################

The equation of the internal energy in non-dimensional form is

.. math::

   \der{T}{t}
   =
   A
   +
   D,

where I introduce some symbols for notational convenience:

.. math::

   A^k
   \equiv
   -
   u_j^k \dder{T^k}{x_j},

.. math::

   D^k
   \equiv
   \frac{1}{\sqrt{Pr} \sqrt{Ra}} \dder{}{x_j} \dder{T^k}{x_j}.

Here the spatial derivatives are approximated by proper finite-difference schemes discussed in :ref:`the other part <spatial_discretisation>`.

**************
Discretisation
**************

The temporal discretisation of the above equation leads to

.. math::

   &\text{do}\,\,\,\, k = 0, 2 \\
   &\,\,\,\,
   \Delta T
   =
   \alpha^k \Delta t \left( A^{k  } + D^{k  } \right)
   +
   \beta^k  \Delta t \left( A^{k-1} + D^{k-1} \right), \\
   &\,\,\,\,
   T^{k+1}
   =
   T^k
   +
   \Delta T. \\
   &\text{enddo}

when all diffusive terms are treated explicitly, while

.. math::

   &\text{do}\,\,\,\, k = 0, 2 \\
   &\,\,\,\,\Delta T
   =
   \alpha^k \Delta t A^{k  }
   +
   \beta^k  \Delta t A^{k-1}
   +
   \gamma^k \Delta t D^{k  }, \\
   &\,\,\,\,T^{k+1}
   =
   T^k
   +
   \left(
      1
      -
      \frac{\gamma^k \Delta t}{2} \frac{1}{\sqrt{Pr} \sqrt{Ra}} \dder{}{y} \dder{}{y}
   \right)^{-1}
   \left(
      1
      -
      \frac{\gamma^k \Delta t}{2} \frac{1}{\sqrt{Pr} \sqrt{Ra}} \dder{}{x} \dder{}{x}
   \right)^{-1}
   \Delta T, \\
   &\text{enddo}

when all diffusive terms are treated implicitly.

Diffusive terms are sometimes partially treated implicitly.
For example, when only the diffusive term in :math:`x` direction is implicitly treated, I have

.. math::

   &\text{do}\,\,\,\, k = 0, 2 \\
   &\,\,\,\,\Delta T
   =
   \alpha^k \Delta t \left( A^{k  } + D_y^{k  } \right)
   +
   \beta^k  \Delta t \left( A^{k-1} + D_y^{k-1} \right)
   +
   \gamma^k \Delta t D_x^{k  }, \\
   &\,\,\,\,T^{k+1}
   =
   T^k
   +
   \left(
      1 - \frac{\gamma^k \Delta t}{2} \frac{1}{\sqrt{Pr} \sqrt{Ra}} \dder{}{x} \dder{}{x}
   \right)^{-1}
   \Delta T, \\
   &\text{enddo}

where

.. math::

   D_x^k
   \equiv
   \frac{1}{\sqrt{Pr} \sqrt{Ra}} \dder{}{x} \dder{T^k}{x},

.. math::

   D_y^k
   \equiv
   \frac{1}{\sqrt{Pr} \sqrt{Ra}} \dder{}{y} \dder{T^k}{y}.

.. seealso::

   :ref:`Time-marching schemes <time_marchers>`.

**************
Implementation
**************

#. Compute right-hand-side terms

   :ref:`fluid_compute_rhs <fluid_compute_rhs>`

#. Compute :math:`\Delta T`

   :ref:`fluid_predict_field <fluid_predict_field>`

#. Solve linear systems and update temperature field

   :ref:`fluid_predict_field <fluid_predict_field>`

