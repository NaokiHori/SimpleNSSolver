
.. _impl_dif_y_z:

############################
Diffusion of y momentum in z
############################

.. math::

   \vat{
      \frac{\sqrt{Pr}}{\sqrt{Ra}} \dder{}{z} \dder{\uy}{z}
   }{\xic, \xjc, \xkc},

which is the product of the diffusivity and

.. math::

   l \vat{\uy}{\xic, \xjc, \xkmm}
   +
   c \vat{\uy}{\xic, \xjc, \xkc }
   +
   u \vat{\uy}{\xic, \xjc, \xkpp},

where

.. math::

   l
   \equiv
   \frac{1}{\Delta z}
   \frac{1}{\Delta z}

.. math::

   u
   \equiv
   \frac{1}{\Delta z}
   \frac{1}{\Delta z}

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
   :tag: uy is diffused in z

.. myliteralinclude:: /../../src/fluid/integrate/uy.c
   :language: c
   :tag: Laplacian in z

