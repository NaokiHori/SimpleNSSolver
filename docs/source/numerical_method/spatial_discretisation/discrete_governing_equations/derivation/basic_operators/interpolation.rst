
.. _basic_operators_interpolation:

#############
Interpolation
#############

Because of :ref:`the staggered grid arrangement <domain_setup>`, variables are sometimes not defined at the locations I want.
In such cases, I need to interpolate the variables, which is discussed here.

Three different symbols are used in this project, whose superscripts indicate the direction of the interpolation.

*****************************************
Arithmetic average :math:`\dintrpa{q}{x}`
*****************************************

This is the most intuitive averaging.

=============================
From cell-center to cell-face
=============================

.. math::

   2 \vat{\dintrpa{q}{x}}{\xic}
   \equiv
   \vat{q}{\xim}
   +
   \vat{q}{\xip}.

=============================
From cell-face to cell-center
=============================

.. math::

   2 \vat{\dintrpa{q}{x}}{\pic}
   \equiv
   \vat{q}{\pim}
   +
   \vat{q}{\pip}.

*************************************
Volume average :math:`\dintrpv{q}{x}`
*************************************

This is the average using the cell sizes as weights.

=============================
From cell-center to cell-face
=============================

.. math::

   2 \Delta x_{\xic} \vat{\dintrpv{q}{x}}{\xic}
   \equiv
   \Delta x_{\xim} \vat{q}{\xim}
   +
   \Delta x_{\xip} \vat{q}{\xip}.

=============================
From cell-face to cell-center
=============================

.. math::

   2 \Delta x_{\pic} \vat{\dintrpv{q}{x}}{\pic}
   \equiv
   \Delta x_{\pim} \vat{q}{\pim}
   +
   \Delta x_{\pip} \vat{q}{\pip}.

.. note::

   :math:`\dintrpv{q}{}` is identical to :math:`\dintrpa{q}{}` if it is used in the homegeneous (:math:`y` or :math:`z`) directions, since I assume that the grid points are equidistantly spaced.

**********************************
Placeholder :math:`\dintrpu{q}{x}`
**********************************

This symbol indicates that it does not have an explicit formulation yet, i.e. I do not know anything except it is interpolated in the given direction (prefix ``u`` implies unknown).
This is used as a placeholder, which must be replaced by one of the above two averages.

