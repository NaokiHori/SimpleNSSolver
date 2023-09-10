
.. _impl_adv_y_z:

############################
Advection of y momentum in z
############################

.. math::

   -
   \vat{
      \dintrpa{
         \dintrpa{\uz}{y}
         \dder{\uy}{z}
      }{z}
   }{\yic, \yjc, \ykc}
   =
   &
   -
   \frac{1}{2}
   \vat{
      \dintrpa{\uz}{y}
   }{\yic, \yjc, \ykm}
   \vat{
      \dder{\uy}{z}
   }{\yic, \yjc, \ykm} \\
   &
   -
   \frac{1}{2}
   \vat{
      \dintrpa{\uz}{y}
   }{\yic, \yjc, \ykp}
   \vat{
      \dder{\uy}{z}
   }{\yic, \yjc, \ykp},

and thus

.. math::

   l_{\yic, \yjc, \ykc} {\uy}_{\yic, \yjc, \ykmm}
   +
   c_{\yic, \yjc, \ykc} {\uy}_{\yic, \yjc, \ykc }
   +
   u_{\yic, \yjc, \ykc} {\uy}_{\yic, \yjc, \ykpp},

where

.. math::

   l_{\yic, \yjc, \ykc}
   \equiv
   +
   \frac{1}{2}
   \frac{1}{\Delta y}
   \vat{
      \dintrpa{\uz}{y}
   }{\yic, \yjc, \ykm},

.. math::

   u_{\yic, \yjc, \ykc}
   \equiv
   -
   \frac{1}{2}
   \frac{1}{\Delta z}
   \vat{
      \dintrpa{\uz}{y}
   }{\yic, \yjc, \ykp},

.. math::

   c_{\yic, \yjc, \ykc}
   \equiv
   -
   l_{\yic, \yjc, \ykc}
   -
   u_{\yic, \yjc, \ykc}.

.. myliteralinclude:: /../../src/fluid/integrate/uy.c
   :language: c
   :tag: uy is transported by uz

