
.. _impl_dif_x_y:

############################
Diffusion of x momentum in y
############################

.. math::

   \vat{
      \frac{\sqrt{Pr}}{\sqrt{Ra}} \dder{}{y} \dder{\ux}{y}
   }{\xic, \xjc, \xkc},

which is the product of the diffusivity and

.. math::

   l {\ux}_{\xic, \xjmm, \xkc}
   +
   c {\ux}_{\xic, \xjc , \xkc}
   +
   u {\ux}_{\xic, \xjpp, \xkc},

where

.. math::

   l
   \equiv
   \frac{1}{\Delta y} \frac{1}{\Delta y},

.. math::

   u
   \equiv
   \frac{1}{\Delta y} \frac{1}{\Delta y},

.. math::

   c
   \equiv
   -
   l
   -
   u.

.. myliteralinclude:: /../../src/fluid/integrate/ux.c
   :language: c
   :tag: ux is diffused in y

.. myliteralinclude:: /../../src/fluid/integrate/ux.c
   :language: c
   :tag: Laplacian in y

