
.. _impl_dif_t_x:

#############################
Diffusion of temperature in x
#############################

.. math::

   \vat{
      \frac{1}{\sqrt{Pr} \sqrt{Ra}} \dder{}{x} \dder{T}{x}
   }{\pic, \pjc, \pkc},

which is the product of the diffusivity and

.. math::

   l_{\pic} T_{\pimm, \pjc, \pkc}
   +
   c_{\pic} T_{\pic , \pjc, \pkc}
   +
   u_{\pic} T_{\pipp, \pjc, \pkc},

where

.. math::

   l_{\pic}
   \equiv
   \frac{1}{\Delta x_{\pic}} \frac{1}{\Delta x_{\pim}},

.. math::

   u_{\pic}
   \equiv
   \frac{1}{\Delta x_{\pic}} \frac{1}{\Delta x_{\pip}},

.. math::

   c_{\pic}
   \equiv
   -
   l_{\pic}
   -
   u_{\pic}.

.. myliteralinclude:: /../../src/fluid/integrate/t.c
   :language: c
   :tag: T is diffused in x

.. myliteralinclude:: /../../src/fluid/integrate/t.c
   :language: c
   :tag: Laplacian w.r.t. temp in x

