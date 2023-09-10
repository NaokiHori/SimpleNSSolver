
.. _impl_adv_x_y:

############################
Advection of x momentum in y
############################

.. math::

   -
   \vat{
      \dintrpa{
         \dintrpv{\uy}{x}
         \dder{\ux}{y}
      }{y}
   }{\xic, \xjc, \xkc}
   =
   &
   -
   \frac{1}{2} \vat{\dintrpv{\uy}{x}}{\xic, \xjm, \xkc} \vat{\dder{\ux}{y}}{\xic, \xjm, \xkc} \\
   &
   -
   \frac{1}{2} \vat{\dintrpv{\uy}{x}}{\xic, \xjp, \xkc} \vat{\dder{\ux}{y}}{\xic, \xjp, \xkc},

and thus

.. math::

   l_{\xic, \xjc, \xkc} {\ux}_{\xic, \xjmm, \xkc}
   +
   c_{\xic, \xjc, \xkc} {\ux}_{\xic, \xjc , \xkc}
   +
   u_{\xic, \xjc, \xkc} {\ux}_{\xic, \xjpp, \xkc},

where

.. math::

   l_{\xic, \xjc, \xkc}
   \equiv
   + \frac{1}{2} \frac{1}{\Delta y} \vat{\dintrpv{\uy}{x}}{\xic, \xjm, \xkc},

.. math::

   u_{\xic, \xjc, \xkc}
   \equiv
   - \frac{1}{2} \frac{1}{\Delta y} \vat{\dintrpv{\uy}{x}}{\xic, \xjp, \xkc},

.. math::

   c_{\xic, \xjc, \xkc}
   \equiv
   -
   l_{\xic, \xjc, \xkc}
   -
   u_{\xic, \xjc, \xkc}.

.. myliteralinclude:: /../../src/fluid/integrate/ux.c
   :language: c
   :tag: ux is transported by uy

