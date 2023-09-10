
.. _impl_adv_y_y:

############################
Advection of y momentum in y
############################

.. math::

   -
   \vat{
      \dintrpa{
         \dintrpa{\uy}{y}
         \dder{\uy}{y}
      }{y}
   }{\yic, \yjc, \ykc}
   =
   &
   -
   \frac{1}{2}
   \vat{
      \dintrpa{\uy}{y}
   }{\yic, \yjm, \ykc}
   \vat{
      \dder{\uy}{y}
   }{\yic, \yjm, \ykc} \\
   &
   -
   \frac{1}{2}
   \vat{
      \dintrpa{\uy}{y}
   }{\yic, \yjp, \ykc}
   \vat{
      \dder{\uy}{y}
   }{\yic, \yjp, \ykc},

and thus

.. math::

   l_{\yic, \yjc, \ykc} {\uy}_{\yic, \yjmm, \ykc}
   +
   c_{\yic, \yjc, \ykc} {\uy}_{\yic, \yjc , \ykc}
   +
   u_{\yic, \yjc, \ykc} {\uy}_{\yic, \yjpp, \ykc},

where

.. math::

   l_{\yic, \yjc, \ykc}
   \equiv
   +
   \frac{1}{2}
   \frac{1}{\Delta y}
   \vat{
      \dintrpa{\uy}{y}
   }{\yic, \yjm, \ykc},

.. math::

   u_{\yic, \yjc, \ykc}
   \equiv
   -
   \frac{1}{2}
   \frac{1}{\Delta y}
   \vat{
      \dintrpa{\uy}{y}
   }{\yic, \yjp, \ykc},

.. math::

   c_{\yic, \yjc, \ykc}
   \equiv
   -
   l_{\yic, \yjc, \ykc}
   -
   u_{\yic, \yjc, \ykc}.

.. myliteralinclude:: /../../src/fluid/integrate/uy.c
   :language: c
   :tag: uy is transported by uy

