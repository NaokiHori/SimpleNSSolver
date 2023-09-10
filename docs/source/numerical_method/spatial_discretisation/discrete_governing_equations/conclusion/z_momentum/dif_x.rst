
.. _impl_dif_z_x:

############################
Diffusion of z momentum in x
############################

.. math::

   \vat{
      \frac{\sqrt{Pr}}{\sqrt{Ra}} \dder{}{x} \dder{\uz}{x}
   }{\zic, \zjc, \zkc},

which is the product of the diffusivity and

.. math::

   l_{\zic} \vat{\uz}{\zimm, \zjc, \zkc}
   +
   c_{\zic} \vat{\uz}{\zic , \zjc, \zkc}
   +
   u_{\zic} \vat{\uz}{\zipp, \zjc, \zkc},

where

.. math::

   \vat{l}{\zic}
   \equiv
   \frac{1}{\vat{\Delta x}{\zic}}
   \frac{1}{\vat{\Delta x}{\zim}},

.. math::

   \vat{u}{\zic}
   \equiv
   \frac{1}{\vat{\Delta x}{\zic}}
   \frac{1}{\vat{\Delta x}{\zip}},

.. math::

   \vat{c}{\zic}
   \equiv
   -
   \vat{l}{\zic}
   -
   \vat{u}{\zic}.

.. myliteralinclude:: /../../src/fluid/integrate/uz.c
   :language: c
   :tag: uz is diffused in x

.. myliteralinclude:: /../../src/fluid/integrate/uz.c
   :language: c
   :tag: Laplacian w.r.t. uz in x

