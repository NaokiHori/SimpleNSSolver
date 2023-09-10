
.. _impl_dif_y_x:

############################
Diffusion of y momentum in x
############################

.. math::

   \vat{
      \frac{\sqrt{Pr}}{\sqrt{Ra}} \dder{}{x} \dder{\uy}{x}
   }{\yic, \yjc, \ykc},

which is the product of the diffusivity and

.. math::

   l_{\yic} \vat{\uy}{\yimm, \yjc, \ykc}
   +
   c_{\yic} \vat{\uy}{\yic , \yjc, \ykc}
   +
   u_{\yic} \vat{\uy}{\yipp, \yjc, \ykc},

where

.. math::

   \vat{l}{\yic}
   \equiv
   \frac{1}{\vat{\Delta x}{\yic}} \frac{1}{\vat{\Delta x}{\yim}},

.. math::

   \vat{u}{\yic}
   \equiv
   \frac{1}{\vat{\Delta x}{\yic}} \frac{1}{\vat{\Delta x}{\yip}},

.. math::

   \vat{c}{\yic}
   \equiv
   -
   \vat{l}{\yic}
   -
   \vat{u}{\yic}.

.. myliteralinclude:: /../../src/fluid/integrate/uy.c
   :language: c
   :tag: uy is diffused in x

.. myliteralinclude:: /../../src/fluid/integrate/uy.c
   :language: c
   :tag: Laplacian w.r.t. uy in x

