
.. _impl_adv_z_x:

############################
Advection of z momentum in x
############################

.. math::

   -
   \vat{
      \dintrpv{
         \dintrpa{\ux}{z}
         \dder{\uz}{x}
      }{x}
   }{\zic, \zjc, \zkc}
   =
   &
   -
   \frac{1}{2}
   \frac{\Delta x_{\zim}}{\Delta x_{\zic}}
   \vat{
      \dintrpa{\ux}{z}
   }{\zim, \zjc, \zkc}
   \vat{
      \dder{\uz}{x}
   }{\zim, \zjc, \zkc} \\
   &
   -
   \frac{1}{2}
   \frac{\Delta x_{\zip}}{\Delta x_{\zic}}
   \vat{
      \dintrpa{\ux}{z}
   }{\zip, \zjc, \zkc}
   \vat{
      \dder{\uz}{x}
   }{\zip, \zjc, \zkc},

and thus

.. math::

   l_{\zic, \zjc, \zkc} {\uz}_{\zimm, \zjc, \zkc}
   +
   c_{\zic, \zjc, \zkc} {\uz}_{\zic , \zjc, \zkc}
   +
   u_{\zic, \zjc, \zkc} {\uz}_{\zipp, \zjc, \zkc},

where

.. math::

   l_{\zic, \zjc, \zkc}
   \equiv
   +
   \frac{1}{2}
   \frac{\Delta x_{\zim}}{\Delta x_{\zic}}
   \vat{\dintrpa{\ux}{z}}{\zim, \zjc, \zkc}
   \frac{1}{\Delta x_{\zim}}
   =
   +
   \frac{1}{2}
   \frac{1}{\Delta x_{\zic}}
   \vat{\dintrpa{\ux}{z}}{\zim, \zjc, \zkc},

.. math::

   u_{\zic, \zjc, \zkc}
   \equiv
   -
   \frac{1}{2}
   \frac{\Delta x_{\zip}}{\Delta x_{\zic}}
   \vat{\dintrpa{\ux}{z}}{\zip, \zjc, \zkc}
   \frac{1}{\Delta x_{\zip}}
   =
   -
   \frac{1}{2}
   \frac{1}{\Delta x_{\zic}}
   \vat{\dintrpa{\ux}{z}}{\zip, \zjc, \zkc},

.. math::

   c_{\zic, \zjc, \zkc}
   \equiv
   -
   l_{\zic, \zjc, \zkc}
   -
   u_{\zic, \zjc, \zkc}.

.. myliteralinclude:: /../../src/fluid/integrate/uz.c
   :language: c
   :tag: uz is transported by ux

