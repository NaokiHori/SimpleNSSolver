In this part, I consider to discretise

.. math::

   \der{\ux}{t}
   +
   u_j \der{\ux}{x_j}
   =
   -
   \der{p}{x}
   +
   \frac{\sqrt{Pr}}{\sqrt{Ra}} \der{}{x_j} \der{\ux}{x_j}
   +
   T

at the :math:`x` cell faces :math:`\left( \xic, \xjc, \xkc \right)`, and

.. math::

   \der{\uy}{t}
   +
   u_j \der{\uy}{x_j}
   =
   -
   \der{p}{y}
   +
   \frac{\sqrt{Pr}}{\sqrt{Ra}} \der{}{x_j} \der{\uy}{x_j}

at the :math:`y` cell faces :math:`\left( \yic, \yjc, \ykc \right)`, and

.. math::

   \der{\uz}{t}
   +
   u_j \der{\uz}{x_j}
   =
   -
   \der{p}{z}
   +
   \frac{\sqrt{Pr}}{\sqrt{Ra}} \der{}{x_j} \der{\uz}{x_j}

at the :math:`z` cell faces :math:`\left( \zic, \zjc, \zkc \right)`.

Also the kinetic energy balance

.. math::

   \der{k}{t}
   +
   u_i \der{k}{x_i}
   =
   -
   u_i \der{p}{x_i}
   +
   \frac{\sqrt{Pr}}{\sqrt{Ra}} \der{}{x_j} \left( u_i \der{u_i}{x_j} \right)
   -
   \frac{\sqrt{Pr}}{\sqrt{Ra}} \der{u_i}{x_j} \der{u_i}{x_j}
   + u_i T \delta_{ix}

is considered, which is dependent on the momentum equation in the continuous domain (see :ref:`the continuous space <governing_equations>`).

