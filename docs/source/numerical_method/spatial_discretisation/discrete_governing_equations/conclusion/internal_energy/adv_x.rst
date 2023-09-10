
.. _impl_adv_t_x:

#############################
Advection of temperature in x
#############################

.. math::

   -
   \vat{
      \dintrpv{
         \ux
         \dder{T}{x}
      }{x}
   }{\pic, \pjc, \pkc}
   =
   &
   -
   \frac{1}{2}
   \frac{\Delta x_{\pim}}{\Delta x_{\pic}}
   \vat{\ux}{\pim, \pjc, \pkc}
   \vat{\dder{T}{x}}{\pim, \pjc, \pkc} \\
   &
   -
   \frac{1}{2}
   \frac{\Delta x_{\pip}}{\Delta x_{\pic}}
   \vat{\ux}{\pip, \pjc, \pkc}
   \vat{\dder{T}{x}}{\pip, \pjc, \pkc},

and thus

.. math::

   l_{\pic, \pjc, \pkc} T_{\pimm, \pjc, \pkc}
   +
   c_{\pic, \pjc, \pkc} T_{\pic , \pjc, \pkc}
   +
   u_{\pic, \pjc, \pkc} T_{\pipp, \pjc, \pkc},

where

.. math::

   l_{\pic, \pjc, \pkc}
   \equiv
   +
   \frac{1}{2}
   \frac{\Delta x_{\pim}}{\Delta x_{\pic}}
   \vat{\ux}{\pim, \pjc, \pkc}
   \frac{1}{\Delta x_{\pim}}
   =
   +
   \frac{1}{2}
   \frac{1}{\Delta x_{\pic}}
   \vat{\ux}{\pim, \pjc, \pkc},

.. math::

   u_{\pic, \pjc, \pkc}
   \equiv
   -
   \frac{1}{2}
   \frac{\Delta x_{\pip}}{\Delta x_{\pic}}
   \vat{\ux}{\pip, \pjc, \pkc}
   \frac{1}{\Delta x_{\pip}}
   =
   -
   \frac{1}{2}
   \frac{1}{\Delta x_{\pic}}
   \vat{\ux}{\pip, \pjc, \pkc},

.. math::

   c_{\pic, \pjc, \pkc}
   \equiv
   -
   l_{\pic, \pjc, \pkc}
   -
   u_{\pic, \pjc, \pkc}.

.. myliteralinclude:: /../../src/fluid/integrate/t.c
   :language: c
   :tag: T is transported by ux

