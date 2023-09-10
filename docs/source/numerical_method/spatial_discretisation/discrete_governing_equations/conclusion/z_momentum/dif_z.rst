
.. _impl_dif_z_z:

############################
Diffusion of z momentum in z
############################

.. math::

   \vat{
      \frac{\sqrt{Pr}}{\sqrt{Ra}} \dder{}{z} \dder{\uz}{z}
   }{\xic, \xjc, \xkc},

which is the product of the diffusivity and

.. math::

   l \vat{\uz}{\xic, \xjc, \xkmm}
   +
   c \vat{\uz}{\xic, \xjc, \xkc }
   +
   u \vat{\uz}{\xic, \xjc, \xkpp},

where

.. math::

   l
   \equiv
   \frac{1}{\Delta z}
   \frac{1}{\Delta z},

.. math::

   u
   \equiv
   \frac{1}{\Delta z}
   \frac{1}{\Delta z},

.. math::

   c
   \equiv
   - l
   - u.

.. myliteralinclude:: /../../src/fluid/integrate/uz.c
   :language: c
   :tag: uz is diffused in z

.. myliteralinclude:: /../../src/fluid/integrate/uz.c
   :language: c
   :tag: Laplacian in z

