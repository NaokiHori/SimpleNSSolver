
.. _impl_dif_x_x:

############################
Diffusion of x momentum in x
############################

.. math::

   \vat{
      \frac{\sqrt{Pr}}{\sqrt{Ra}} \dder{}{x}  \dder{\ux}{x}
   }{\xic, \xjc, \xkc},

which is the product of the diffusivity and

.. math::

   l_{\xic} {\ux}_{\ximm, \xjc, \xkc}
   +
   c_{\xic} {\ux}_{\xic , \xjc, \xkc}
   +
   u_{\xic} {\ux}_{\xipp, \xjc, \xkc},

where

.. math::

   l_{\xic}
   \equiv
   \frac{1}{\Delta x_{\xic}} \frac{1}{\Delta x_{\xim}},

.. math::

   u_{\xic}
   \equiv
   \frac{1}{\Delta x_{\xic}} \frac{1}{\Delta x_{\xip}},

.. math::

   c_{\xic}
   \equiv
   -
   l_{\xic}
   -
   u_{\xic}.

.. myliteralinclude:: /../../src/fluid/integrate/ux.c
   :language: c
   :tag: ux is diffused in x

.. myliteralinclude:: /../../src/fluid/integrate/ux.c
   :language: c
   :tag: Laplacian w.r.t. ux in x

