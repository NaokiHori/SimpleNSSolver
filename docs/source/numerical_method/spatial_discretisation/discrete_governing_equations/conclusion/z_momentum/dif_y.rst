
.. _impl_dif_z_y:

############################
Diffusion of z momentum in y
############################

.. math::

   \vat{
      \frac{\sqrt{Pr}}{\sqrt{Ra}} \dder{}{y} \dder{\uy}{y}
   }{\zic, \zjc, \zkc},

which is the product of the diffusivity and

.. math::

   l \vat{\uy}{\zic, \zjmm, \zkc}
   +
   c \vat{\uy}{\zic, \zjc , \zkc}
   +
   u \vat{\uy}{\zic, \zjpp, \zkc},

where

.. math::

   l
   \equiv
   \frac{1}{\Delta y}
   \frac{1}{\Delta y},

.. math::

   u
   \equiv
   \frac{1}{\Delta y}
   \frac{1}{\Delta y},

.. math::

   c
   \equiv
   - l
   - u.

.. myliteralinclude:: /../../src/fluid/integrate/uz.c
   :language: c
   :tag: uz is diffused in y

.. myliteralinclude:: /../../src/fluid/integrate/uz.c
   :language: c
   :tag: Laplacian in y

