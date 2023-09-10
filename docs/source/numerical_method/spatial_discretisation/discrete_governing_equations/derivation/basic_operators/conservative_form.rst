
.. _basic_operators_conservative:

#################
Conservative form
#################

To begin with, I provide a definition of the terminology *conservative form* in both the continuous and discrete domains.
For simplicity, I assume that the domain is one-dimensional

.. math::

   x \in \left[ 0, l_x \right],

whose edges (:math:`x = 0`, :math:`x = l_x`) are bounded or periodic.

In this simple set-up, i consider

.. math::

   \der{q}{x}.

Integrating this term in the whole domain

.. math::

   \int_{x = 0}^{x = l_x} \der{q}{x} dx

of course yields

.. math::

   \vat{q}{l_x} - \vat{q}{0}.

This equation states that the change in :math:`q` is determined by the fluxes at the boundaries, and if the fluxes are zero (or the domain is periodic), the total amount of :math:`q` is *conserved*.

Thus, in the continuous domain, I call the term is *conservative* (or in the *conservative form*) when it is written as

.. math::

   \der{q}{x_i},

where :math:`q` can be an arbitrary tensor of any order.

