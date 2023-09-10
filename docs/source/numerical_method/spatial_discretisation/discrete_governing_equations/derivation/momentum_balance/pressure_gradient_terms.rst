#######################
Pressure-gradient terms
#######################

********
Momentum
********

Since they are the simple gradients, I discretise them as

.. math::

   \vat{\dder{p}{x}}{\xic, \xjc},
   \vat{\dder{p}{y}}{\yic, \yjc}

in both directions, respectively.

****************
Squared velocity
****************

I have

.. math::

   \vat{\ux \dder{p}{x}}{\xic, \xjc}

and

.. math::

   \vat{\uy \dder{p}{y}}{\yic, \yjc}

in each direction.

To see the global effect on the squared velocity balances, I consider to discretely integrate these terms in the whole volume.
Fist I focus on the :math:`x` contribution:

.. math::

   \sum_{\xjc} \sum_{\xic}
   \vat{\ux \dder{p}{x}}{\xic, \xjc}
   \Delta x_{\xic}
   \Delta y_{\xjc}
   =
   \sum_{\xjc} \sum_{\xic}
   \vat{\ux}{\xic, \xjc}
   \left(
      \vat{p}{\xip, \xjc}
      -
      \vat{p}{\xim, \xjc}
   \right)
   \Delta y_{\xjc}.

Although this equation is written as the sum of each :math:`\ux`, I can write it as the sum of each :math:`p`:

.. math::

   -
   \sum_{\pjc} \sum_{\pic}
   \vat{p}{\pic, \pjc}
   \left(
      \vat{\ux}{\pip, \pjc}
      -
      \vat{\ux}{\pim, \pjc}
   \right)
   \Delta y_{\pjc},

where I use the impermeable condition to eliminate the boundary contributions.

Similarly, the :math:`y` contribution yields

.. math::

   -
   \sum_{\pjc} \sum_{\pic}
   \vat{p}{\pic, \pjc}
   \left(
      \vat{\uy}{\pic, \pjp}
      -
      \vat{\uy}{\pic, \pjm}
   \right)
   \Delta x_{\pjc}.

Although they do not vanish independently, when added:

.. math::

   &
   -
   \sum_{\pjc} \sum_{\pic}
   \vat{p}{\pic, \pjc}
   \left(
      \vat{\ux}{\pip, \pjc}
      -
      \vat{\ux}{\pim, \pjc}
   \right)
   \Delta y_{\pjc}
   -
   \sum_{\pjc} \sum_{\pic}
   \vat{p}{\pic, \pjc}
   \left(
      \vat{\uy}{\pic, \pjp}
      -
      \vat{\uy}{\pic, \pjm}
   \right)
   \Delta x_{\pjc} \\
   =
   &
   -
   \sum_{\pjc} \sum_{\pic}
   \vat{p}{\pic, \pjc}
   \left[
      \frac{
         \vat{\ux}{\pip, \pjc}
         -
         \vat{\ux}{\pim, \pjc}
      }{\Delta x_{\pic}}
      +
      \frac{
         \vat{\uy}{\pic, \pjp}
         -
         \vat{\uy}{\pic, \pjm}
      }{\Delta y_{\pjc}}
   \right]
   \Delta x_{\pic} \Delta y_{\pjc},

they disappear because of the discrete incompressibility constraint.

