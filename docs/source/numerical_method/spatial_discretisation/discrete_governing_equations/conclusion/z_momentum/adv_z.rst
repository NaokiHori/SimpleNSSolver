
.. _impl_adv_z_z:

############################
Advection of z momentum in z
############################

.. math::

   -
   \vat{
      \dintrpa{
         \dintrpa{\uz}{z}
         \dder{\uz}{z}
      }{z}
   }{\zic, \zjc, \zkc}
   =
   &
   -
   \vat{
      \dintrpa{\uz}{z}
   }{\zic, \zjc, \zkm}
   \vat{
      \dder{\uz}{z}
   }{\zic, \zjc, \zkm} \\
   &
   -
   \vat{
      \dintrpa{\uz}{z}
   }{\zic, \zjc, \zkp}
   \vat{
      \dder{\uz}{z}
   }{\zic, \zjc, \zkp},

and thus

.. math::

   l_{\zic, \zjc, \zkc} {\uz}_{\zic, \zjc, \zkmm}
   +
   c_{\zic, \zjc, \zkc} {\uz}_{\zic, \zjc, \zkc }
   +
   u_{\zic, \zjc, \zkc} {\uz}_{\zic, \zjc, \zkpp},

where

.. math::

   l_{\zic, \zjc, \zkc}
   \equiv
   +
   \frac{1}{2}
   \frac{1}{\Delta z}
   \vat{\dintrpa{\uz}{z}}{\zic, \zjc, \zkm},

.. math::

   u_{\zic, \zjc, \zkc}
   \equiv
   -
   \frac{1}{2}
   \frac{1}{\Delta z}
   \vat{\dintrpa{\uz}{z}}{\zic, \zjc, \zkp},

.. math::

   c_{\zic, \zjc, \zkc}
   \equiv
   -
   l_{\zic, \zjc, \zkc}
   -
   u_{\zic, \zjc, \zkc}.

.. myliteralinclude:: /../../src/fluid/integrate/uz.c
   :language: c
   :tag: uz is transported by uz

