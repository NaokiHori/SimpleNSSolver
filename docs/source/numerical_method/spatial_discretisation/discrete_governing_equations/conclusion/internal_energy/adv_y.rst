
.. _impl_adv_t_y:

#############################
Advection of temperature in y
#############################

.. math::

   -
   \vat{
      \dintrpa{
         \uy
         \dder{T}{y}
      }{y}
   }{\pic, \pjc, \pkc}
   =
   &
   -
   \frac{1}{2}
   \vat{\uy}{\pic, \pjm, \pkc}
   \vat{\dder{T}{y}}{\pic, \pjm, \pkc} \\
   &
   -
   \frac{1}{2}
   \vat{\uy}{\pic, \pjp, \pkc}
   \vat{\dder{T}{y}}{\pic, \pjp, \pkc},

and thus

.. math::

   l_{\pic, \pjc, \pkc} T_{\pic, \pjmm, \pkc}
   +
   c_{\pic, \pjc, \pkc} T_{\pic, \pjc , \pkc}
   +
   u_{\pic, \pjc, \pkc} T_{\pic, \pjpp, \pkc},

where

.. math::

   l_{\pic, \pjc, \pkc}
   \equiv
   +
   \frac{1}{2}
   \frac{1}{\Delta y}
   \vat{\uy}{\pic, \pjm, \pkc},

.. math::

   u_{\pic, \pjc, \pkc}
   \equiv
   -
   \frac{1}{2}
   \frac{1}{\Delta y}
   \vat{\uy}{\pic, \pjp, \pkc},

.. math::

   c_{\pic, \pjc, \pkc}
   \equiv
   -
   l_{\pic, \pjc, \pkc}
   -
   u_{\pic, \pjc, \pkc}.

.. myliteralinclude:: /../../src/fluid/integrate/t.c
   :language: c
   :tag: T is transported by uy

