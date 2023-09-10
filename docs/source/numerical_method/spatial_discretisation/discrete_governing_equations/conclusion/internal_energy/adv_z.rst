
.. _impl_adv_t_z:

#############################
Advection of temperature in z
#############################

.. math::

   -
   \vat{
      \dintrpa{
         \uz
         \dder{T}{z}
      }{z}
   }{\pic, \pjc, \pkc}
   =
   &
   -
   \frac{1}{2}
   \vat{\uz}{\pic, \pjc, \pkm}
   \vat{\dder{T}{z}}{\pic, \pjc, \pkm} \\
   &
   -
   \frac{1}{2}
   \vat{\uz}{\pic, \pjc, \pkp}
   \vat{\dder{T}{z}}{\pic, \pjc, \pkp},

and thus

.. math::

   l_{\pic, \pjc, \pkc} T_{\pic, \pjc, \pkmm}
   +
   c_{\pic, \pjc, \pkc} T_{\pic, \pjc, \pkc }
   +
   u_{\pic, \pjc, \pkc} T_{\pic, \pjc, \pkpp},

where

.. math::

   l_{\pic, \pjc, \pkc}
   \equiv
   +
   \frac{1}{2}
   \frac{1}{\Delta z}
   \vat{\uz}{\pic, \pjc, \pkm},

.. math::

   u_{\pic, \pjc, \pkc}
   \equiv
   -
   \frac{1}{2}
   \frac{1}{\Delta z}
   \vat{\uz}{\pic, \pjc, \pkp},

.. math::

   c_{\pic, \pjc, \pkc}
   \equiv
   -
   l_{\pic, \pjc, \pkc}
   -
   u_{\pic, \pjc, \pkc}.

.. myliteralinclude:: /../../src/fluid/integrate/t.c
   :language: c
   :tag: T is transported by uz

