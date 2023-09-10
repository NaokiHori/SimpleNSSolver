
.. _impl_adv_y_x:

############################
Advection of y momentum in x
############################

.. math::

   -
   \vat{
      \dintrpv{
         \dintrpa{\ux}{y}
         \dder{\uy}{x}
      }{x}
   }{\yic, \yjc, \ykc}
   =
   &
   -
   \frac{1}{2}
   \frac{\Delta x_{\yim}}{\Delta x_{\yic}}
   \vat{
      \dintrpa{\ux}{y}
   }{\yim, \yjc, \ykc}
   \vat{
      \dder{\uy}{x}
   }{\yim, \yjc, \ykc} \\
   &
   -
   \frac{1}{2} \frac{\Delta x_{\yip}}{\Delta x_{\yic}}
   \vat{
      \dintrpa{\ux}{y}
   }{\yip, \yjc, \ykc}
   \vat{
      \dder{\uy}{x}
   }{\yip, \yjc, \ykc},

and thus

.. math::

   l_{\yic, \yjc, \ykc} {\uy}_{\yimm, \yjc, \ykc}
   +
   c_{\yic, \yjc, \ykc} {\uy}_{\yic , \yjc, \ykc}
   +
   u_{\yic, \yjc, \ykc} {\uy}_{\yipp, \yjc, \ykc},

where

.. math::

   l_{\yic, \yjc, \ykc}
   \equiv
   +
   \frac{1}{2}
   \frac{\Delta x_{\yim}}{\Delta x_{\yic}}
   \vat{\dintrpa{\ux}{y}}{\yim, \yjc, \ykc}
   \frac{1}{\Delta x_{\yim}}
   =
   +
   \frac{1}{2}
   \frac{1}{\Delta x_{\yic}}
   \vat{\dintrpa{\ux}{y}}{\yim, \yjc, \ykc},

.. math::

   u_{\yic, \yjc, \ykc}
   \equiv
   -
   \frac{1}{2}
   \frac{\Delta x_{\yip}}{\Delta x_{\yic}}
   \vat{\dintrpa{\ux}{y}}{\yip, \yjc, \ykc}
   \frac{1}{\Delta x_{\yip}}
   =
   -
   \frac{1}{2}
   \frac{1}{\Delta x_{\yic}}
   \vat{\dintrpa{\ux}{y}}{\yip, \yjc, \ykc},

.. math::

   c_{\yic, \yjc, \ykc}
   \equiv
   -
   l_{\yic, \yjc, \ykc}
   -
   u_{\yic, \yjc, \ykc}.

.. myliteralinclude:: /../../src/fluid/integrate/uy.c
   :language: c
   :tag: uy is transported by ux

