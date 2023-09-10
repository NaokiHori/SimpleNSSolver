############################
Differentiation and integral
############################

.. note::

   Hereafter :math:`\der{q}{x}` (\\partial) and :math:`\dder{q}{x}` (\\delta) are used to distinguish the derivatives in the continuous and discrete domains, respectively.

To mimic the property of the conservative form in the discrete domain, I should discretise the first-order derivatives as follows.

******************************
Quantity defined at cell faces
******************************

When the quantity is defined at cell faces, its first-order derivative should be defined at cell centers:

.. math::

   \vat{\dder{q}{x}}{\pic}
   =
   \frac{\vat{q}{\pip} - \vat{q}{\pim}}{\vat{x}{\pip} - \vat{x}{\pim}}.

At the same time, the integral of a quantity which is defined at cell centers should be discretised as

.. math::

   \sum_{center}
   \vat{q}{\pic}
   \vat{\Delta x}{\pic},

or hereafter simply

.. math::

   \sum_{\pic}
   q
   \Delta x,

where cell-centered grid sizes are used to approximate :math:`dx`.

.. mydetails:: Derivation

   .. note::

      Although the derivation below is so verbase, I show it here for completeness.

   I consider a quantity :math:`q`, whose derivative is defined at each cell center, i.e.

   .. math::

      \vat{\dder{q}{x}}{\pic}.

   Integrating the term in the whole domain discretely yields

   .. math::

      \int \der{q}{x} dx
      \approx
      \sum_{i} \vat{
         \left(
            \dder{q}{x} \Delta x
         \right)
      }{\pic}.

   Since :math:`\vat{\Delta x}{\pic}` denotes the cell size defined at each cell center, it is clearly

   .. math::

      \vat{\Delta x}{\pic}
      =
      \vat{x}{\pip}
      -
      \vat{x}{\pim}.

   Now I consider to *mimic* the relation which is satisfied in the continuous domain: namely, the volume integral of this term should be given by the boundary values.
   To do so, obviously the denominator of the gradient :math:`\delta x` should be equal to the size of the cell :math:`\vat{\Delta x}{\pic}`, giving

   .. math::

      \vat{\delta x}{\pic}
      =
      \vat{x}{\pip}
      -
      \vat{x}{\pim}.

   Also, it is natural to use the same scheme for the numerator :math:`\delta q` since they have the same form:

   .. math::

      \vat{\delta q}{\pic}
      =
      \vat{q}{\pip}
      -
      \vat{q}{\pim}.

   Thus, I conclude that

   .. math::

      \vat{\dder{q}{x}}{\pic}

   should be described as

   .. math::

      \frac{\vat{q}{\pip} - \vat{q}{\pim}}{\vat{x}{\pip} - \vat{x}{\pim}}.

********************************
Quantity defined at cell centers
********************************

When the quantity is defined at cell centers (in general, see below), its first-order derivative should be defined at cell faces:

.. math::

   \vat{\dder{q}{x}}{\pip}
   =
   \frac{\vat{q}{\pipp} - \vat{q}{\pic}}{\vat{x}{\pipp} - \vat{x}{\pic}}.

At the same time, the integral of a quantity which is defined at cell faces should be discretised as

.. math::

   \sum_{face}
   \vat{q}{\pip}
   \vat{\Delta x}{\pip},

or hereafter simply

.. math::

   \sum_{\pip}
   q
   \Delta x,

where cell-faced grid sizes are used to approximate :math:`dx`.

The derivation is omitted since it is similar to the one for the cell-centered quantity discussed above.

.. note::

   In the vicinity of the boundaries, the cell-centered quantities are defined on the boundaries, i.e. at the cell-face positions.
   Please refer to :ref:`the domain set-up <domain_setup>`.

***********************************************
Inconsistent (non-conservative) differentiation
***********************************************

Sometimes the first-order derivative at cell center is given by

.. math::

   \frac{\vat{q}{\pipp} - \vat{q}{\pimm}}{\vat{x}{\pipp} - \vat{x}{\pimm}},

i.e. using the neighbouring cell-center values to evaluate the differentiation at a cell center.

To keep the relation between the differentiation and the integral, the volume integral should be approximated correspondingly:

.. math::

   \int q dx
   \approx
   \sum q \left( \vat{x}{\pipp} - \vat{x}{\pimm} \right).

This does not guarantee the conservative nature, since the neighbouring :math:`q` (e.g. :math:`q_{\pic}` and :math:`q_{\pipp}`) are not related and thus the bulk values do not cancel out to each other.
I should avoid using this kind of inconsistent schemes.

