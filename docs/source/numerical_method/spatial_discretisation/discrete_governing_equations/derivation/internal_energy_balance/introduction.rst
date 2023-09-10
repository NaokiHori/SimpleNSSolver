In this part, I consider to discretise

.. math::

   \der{T}{t}
   +
   u_i \der{T}{x_i}
   =
   \frac{1}{\sqrt{Pr} \sqrt{Ra}} \der{}{x_i} \der{T}{x_i}

at the cell center :math:`\left( \pic, \pjc, \pkc \right)`.

Also the balance of the quadratic quantity :math:`h \equiv T^2 / 2`:

.. math::

   \der{h}{t}
   +
   u_i \der{h}{x_i}
   =
   \frac{1}{\sqrt{Pr} \sqrt{Ra}} \der{}{x_i} \left( T \der{T}{x_i} \right)
   -
   \frac{1}{\sqrt{Pr} \sqrt{Ra}} \der{T}{x_i} \der{T}{x_i}

is considered.

.. note::

   In the continuous domain, I obtain this equation by multiplying the temperature :math:`T` and the equation of the thermal energy balance, which is derived in :ref:`the governing equations <governing_equations>`.

