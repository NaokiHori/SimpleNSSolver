
.. _impl_adv_x_z:

############################
Advection of x momentum in z
############################

.. math::

   -
   \vat{
      \dintrpa{
         \dintrpv{\uz}{x}
         \dder{\ux}{z}
      }{z}
   }{\xic, \xjc, \xkc}
   =
   &
   -
   \frac{1}{2} \vat{\dintrpv{\uz}{x}}{\xic, \xjc, \xkm} \vat{\dder{\ux}{z}}{\xic, \xjc, \xkm} \\
   &
   -
   \frac{1}{2} \vat{\dintrpv{\uz}{x}}{\xic, \xjc, \xkp} \vat{\dder{\ux}{z}}{\xic, \xjc, \xkp},

and thus

.. math::

   l_{\xic, \xjc, \xkc} {\ux}_{\xic, \xjc, \xkmm}
   +
   c_{\xic, \xjc, \xkc} {\ux}_{\xic, \xjc, \xkc }
   +
   u_{\xic, \xjc, \xkc} {\ux}_{\xic, \xjc, \xkpp},

where

.. math::

   l_{\xic, \xjc, \xkc}
   \equiv
   + \frac{1}{2} \frac{1}{\Delta z} \vat{\dintrpv{\uz}{x}}{\xic, \xjc, \xkm},

.. math::

   u_{\xic, \xjc, \xkc}
   \equiv
   - \frac{1}{2} \frac{1}{\Delta z} \vat{\dintrpv{\uz}{x}}{\xic, \xjc, \xkp},

.. math::

   c_{\xic, \xjc, \xkc}
   \equiv
   -
   l_{\xic, \xjc, \xkc}
   -
   u_{\xic, \xjc, \xkc}.

.. myliteralinclude:: /../../src/fluid/integrate/ux.c
   :language: c
   :tag: ux is transported by uz

