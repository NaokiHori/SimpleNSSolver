
.. _impl_adv_z_y:

############################
Advection of z momentum in y
############################

.. math::

   -
   \vat{
      \dintrpa{
         \dintrpa{\uy}{z}
         \dder{\uz}{y}
      }{y}
   }{\zic, \zjc, \zkc}
   =
   &
   -
   \vat{
      \dintrpa{\uy}{z}
   }{\zic, \zjm, \zkc}
   \vat{
      \dder{\uz}{y}
   }{\zic, \zjm, \zkc} \\
   &
   -
   \vat{
      \dintrpa{\uy}{z}
   }{\zic, \zjp, \zkc}
   \vat{
      \dder{\uz}{y}
   }{\zic, \zjp, \zkc},

and thus

.. math::

   l_{\zic, \zjc, \zkc} {\uz}_{\zic, \zjmm, \zkc}
   +
   c_{\zic, \zjc, \zkc} {\uz}_{\zic, \zjc , \zkc}
   +
   u_{\zic, \zjc, \zkc} {\uz}_{\zic, \zjpp, \zkc},

where

.. math::

   l_{\zic, \zjc, \zkc}
   \equiv
   +
   \frac{1}{2}
   \frac{1}{\Delta y}
   \vat{\dintrpa{\uy}{z}}{\zic, \zjm, \zkc},

.. math::

   u_{\zic, \zjc, \zkc}
   \equiv
   -
   \frac{1}{2}
   \frac{1}{\Delta y}
   \vat{\dintrpa{\uy}{z}}{\zic, \zjp, \zkc},

.. math::

   c_{\zic, \zjc, \zkc}
   \equiv
   -
   l_{\zic, \zjc, \zkc}
   -
   u_{\zic, \zjc, \zkc}.

.. myliteralinclude:: /../../src/fluid/integrate/uz.c
   :language: c
   :tag: uz is transported by uy

