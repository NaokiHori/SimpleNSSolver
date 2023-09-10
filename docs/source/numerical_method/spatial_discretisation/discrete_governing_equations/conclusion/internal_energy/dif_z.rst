
.. _impl_dif_t_z:

#############################
Diffusion of temperature in z
#############################

.. math::

   \vat{
      \frac{1}{\sqrt{Pr} \sqrt{Ra}} \dder{}{z} \dder{T}{z}
   }{\pic, \pjc, \pkc},

which is the product of the diffusivity and

.. math::

   l T_{\pic, \pjc, \pkmm}
   +
   c T_{\pic, \pjc, \pkc }
   +
   u T_{\pic, \pjc, \pkpp},

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

.. myliteralinclude:: /../../src/fluid/integrate/t.c
   :language: c
   :tag: T is diffused in z

.. myliteralinclude:: /../../src/fluid/integrate/t.c
   :language: c
   :tag: Laplacian in z

