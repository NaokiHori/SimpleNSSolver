
.. _impl_adv_x_x:

############################
Advection of x momentum in x
############################

.. math::

   -
   \vat{
      \dintrpv{
         \dintrpa{\ux}{x}
         \dder{\ux}{x}
      }{x}
   }{\xic, \xjc, \xkc}
   =
   &
   -
   \frac{1}{2}
   \frac{\Delta x_{\xim}}{\Delta x_{\xic}}
   \vat{\dintrpa{\ux}{x}}{\xim, \xjc, \xkc}
   \vat{\dder{\ux}{x}}{\xim, \xjc, \xkc} \\
   &
   -
   \frac{1}{2}
   \frac{\Delta x_{\xip}}{\Delta x_{\xic}}
   \vat{\dintrpa{\ux}{x}}{\xip, \xjc, \xkc}
   \vat{\dder{\ux}{x}}{\xip, \xjc, \xkc},

and thus

.. math::

   l_{\xic, \xjc, \xkc} {\ux}_{\ximm, \xjc, \xkc}
   +
   c_{\xic, \xjc, \xkc} {\ux}_{\xic , \xjc, \xkc}
   +
   u_{\xic, \xjc, \xkc} {\ux}_{\xipp, \xjc, \xkc},

where

.. math::

   l_{\xic, \xjc, \xkc}
   \equiv
   + \frac{1}{2} \frac{\Delta x_{\xim}}{\Delta x_{\xic}} \vat{\dintrpa{\ux}{x}}{\xim, \xjc, \xkc} \frac{1}{\Delta x_{\xim}}
   =
   + \frac{1}{2} \frac{1}{\Delta x_{\xic}} \vat{\dintrpa{\ux}{x}}{\xim, \xjc, \xkc},

.. math::

   u_{\xic, \xjc, \xkc}
   \equiv
   - \frac{1}{2} \frac{\Delta x_{\xip}}{\Delta x_{\xic}} \vat{\dintrpa{\ux}{x}}{\xip, \xjc, \xkc} \frac{1}{\Delta x_{\xip}}
   =
   - \frac{1}{2} \frac{1}{\Delta x_{\xic}} \vat{\dintrpa{\ux}{x}}{\xip, \xjc, \xkc},

.. math::

   c_{\xic, \xjc, \xkc}
   \equiv
   -
   l_{\xic, \xjc, \xkc}
   -
   u_{\xic, \xjc, \xkc}.

.. myliteralinclude:: /../../src/fluid/integrate/ux.c
   :language: c
   :tag: ux is transported by ux

