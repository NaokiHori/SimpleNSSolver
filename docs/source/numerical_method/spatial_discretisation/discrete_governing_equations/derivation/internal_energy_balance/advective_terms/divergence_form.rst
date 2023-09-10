###############
Divergence form
###############

*******
Summary
*******

I find the discrete form of the advection of :math:`T` is given by

.. math::

   \dder{\ux \dintrpu{T}{x}}{x}
   +
   \dder{\uy \dintrpa{T}{y}}{y},

where :math:`\dintrpu{T}{x}` has not been defined yet.

**********
Derivation
**********

As derived in :ref:`the governing equations <governing_equations>`, the advection of the temperature is described as

.. math::

   u_i \der{T}{x_i},

which is the gradient form of the advective terms.
Thanks to the incompressibility constraint, I can make this term conservative:

.. math::

   u_i \der{T}{x_i}
   +
   T \der{u_i}{x_i}
   =
   \der{u_i T}{x_i},

whose right-hand-side term is known as the divergence form.

Since the divergence form is inherently conservative, I obtain a discrete form of the advective term as

.. math::

   \dder{u_i T}{x_i}
   =
   \dder{\ux \dintrpu{T}{x}}{x}
   +
   \dder{\uy \dintrpu{T}{y}}{y},

which is defined at the cell center :math:`\left( \pic, \pjc \right)`.

.. note::

   Interpolations are not necessary for the velocities :math:`\ux` and :math:`\uy` since they are already there.

Because of the uniform grid spacings in the :math:`y` direction, I can replace the interpolation in the :math:`y` direction with the arithmetic average:

.. math::

   \dder{u_i T}{x_i}
   =
   \dder{\ux \dintrpu{T}{x}}{x}
   +
   \dder{\uy \dintrpa{T}{y}}{y}.

To complete the interpolation in the :math:`x` direction, I need to investigate the advection of the quadratic quantity.

