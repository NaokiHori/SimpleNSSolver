
.. _impl_dif_y_y:

############################
Diffusion of y momentum in y
############################

.. math::

   \vat{
      \frac{\sqrt{Pr}}{\sqrt{Ra}} \dder{}{y} \dder{\uy}{y}
   }{\yic, \yjc, \ykc},

which is the product of the diffusivity and

.. math::

   l \vat{\uy}{\yic, \yjmm, \ykc}
   +
   c \vat{\uy}{\yic, \yjc , \ykc}
   +
   u \vat{\uy}{\yic, \yjpp, \ykc},

where

.. math::

   l
   \equiv
   \frac{1}{\Delta y}
   \frac{1}{\Delta y}

.. math::

   u
   \equiv
   \frac{1}{\Delta y}
   \frac{1}{\Delta y}

and

.. math::

   c
   \equiv
   -
   l
   -
   u.

.. myliteralinclude:: /../../src/fluid/integrate/uy.c
   :language: c
   :tag: uy is diffused in y

.. myliteralinclude:: /../../src/fluid/integrate/uy.c
   :language: c
   :tag: Laplacian in y

