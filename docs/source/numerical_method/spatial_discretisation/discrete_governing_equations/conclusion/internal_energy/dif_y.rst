
.. _impl_dif_t_y:

#############################
Diffusion of temperature in y
#############################

.. math::

   \vat{
      \frac{1}{\sqrt{Pr} \sqrt{Ra}} \dder{}{y} \dder{T}{y}
   }{\pic, \pjc, \pkc},

which is the product of the diffusivity and

.. math::

   l T_{\pic, \pjmm, \pkc}
   +
   c T_{\pic, \pjc , \pkc}
   +
   u T_{\pic, \pjpp, \pkc},

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

.. myliteralinclude:: /../../src/fluid/integrate/t.c
   :language: c
   :tag: T is diffused in y

.. myliteralinclude:: /../../src/fluid/integrate/t.c
   :language: c
   :tag: Laplacian in y

