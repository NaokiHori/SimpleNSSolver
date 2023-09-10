####################
Quadratic quantities
####################

.. mydeclare:: /../../src/logging/energy.c
   :language: c
   :tag: logging_check_energy

This function computes the discrete kinetic energies (quadratic quantities in each direction):

.. math::

   k_x = \int \frac{1}{2} \ux^2 dx dy dz, \\
   k_y = \int \frac{1}{2} \uy^2 dx dy dz, \\
   k_z = \int \frac{1}{2} \uz^2 dx dy dz, \\

and the discrete thermal energy (quadratic quantity):

.. math::

   h = \int \frac{1}{2} T^2 dx dy dz

and outputs the results to a file.

*************************************************
Kinetic quadratic quantity in :math:`x` direction
*************************************************

   .. math::

      \sum_{i} \sum_{j} \sum_{k}
      \vat{
         \left(
            \frac{1}{2}
            \ux^2
            \Delta x
            \Delta y
            \Delta z
         \right)
      }{\xic, \xjc, \xkc}

   .. myliteralinclude:: /../../src/logging/energy.c
      :language: c
      :tag: compute quadratic quantity in x direction

*************************************************
Kinetic quadratic quantity in :math:`y` direction
*************************************************

   .. math::

      \sum_{i} \sum_{j} \sum_{k}
      \vat{
         \left(
            \frac{1}{2}
            \uy^2
            \Delta x
            \Delta y
            \Delta z
         \right)
      }{\yic, \yjc, \ykc}

   .. myliteralinclude:: /../../src/logging/energy.c
      :language: c
      :tag: compute quadratic quantity in y direction

*************************************************
Kinetic quadratic quantity in :math:`z` direction
*************************************************

   .. math::

      \sum_{i} \sum_{j} \sum_{k}
      \vat{
         \left(
            \frac{1}{2}
            \uz^2
            \Delta x
            \Delta y
            \Delta z
         \right)
      }{\zic, \zjc, \zkc}

   .. myliteralinclude:: /../../src/logging/energy.c
      :language: c
      :tag: compute quadratic quantity in z direction

**************************
Thermal quadratic quantity
**************************

   .. math::

      \sum_{i} \sum_{j} \sum_{k}
      \vat{
         \left(
            \frac{1}{2}
            T^2
            \Delta x
            \Delta y
            \Delta z
         \right)
      }{\pic, \pjc, \pkc}

   .. myliteralinclude:: /../../src/logging/energy.c
      :language: c
      :tag: compute thermal energy

