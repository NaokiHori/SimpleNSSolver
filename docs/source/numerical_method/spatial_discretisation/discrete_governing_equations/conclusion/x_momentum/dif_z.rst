
.. _impl_dif_x_z:

############################
Diffusion of x momentum in z
############################

.. math::

   \vat{
      \frac{\sqrt{Pr}}{\sqrt{Ra}} \dder{}{z} \dder{\ux}{z}
   }{\xic, \xjc, \xkc},

which is the product of the diffusivity and

.. math::

   l {\ux}_{\xic, \xjc, \xkmm}
   +
   c {\ux}_{\xic, \xjc, \xkc }
   +
   u {\ux}_{\xic, \xjc, \xkpp},

where

.. math::

   l
   \equiv
   \frac{1}{\Delta z} \frac{1}{\Delta z},

.. math::

   u
   \equiv
   \frac{1}{\Delta z} \frac{1}{\Delta z},

.. math::

   c
   \equiv
   -
   l
   -
   u.

.. myliteralinclude:: /../../src/fluid/integrate/ux.c
   :language: c
   :tag: ux is diffused in z

.. myliteralinclude:: /../../src/fluid/integrate/ux.c
   :language: c
   :tag: Laplacian in z

