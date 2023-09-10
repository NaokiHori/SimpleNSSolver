########################
Derivation: mass balance
########################

.. include:: tilde_note.rst

I consider a closed control volume :math:`V`, whose surface is denoted by :math:`\partial V` or :math:`S` and its outward normal is given by :math:`n_i`.

By assuming there is no mass source nor sink, the mass balance inside this control volume leads to

.. math::

   \der{}{t} \int_{V} \rho dV
   +
   \int_{\partial V} \rho u_i n_i dS
   =
   0,

i.e. the increase or decrease of the mass is determined by the mass flux.
Here, :math:`\rho` and :math:`u_i` denote the density and the velocity of the liquid, which are functions of space :math:`x_i` and time :math:`t`.

Assuming the temporal derivative and the spatial integral are commutative, I have

.. math::

   \int_{V} \der{\rho}{t} dV
   +
   \int_{\partial V} \rho u_i n_i dS
   =
   0,

which yields

.. math::

   \int_{V} \left( \der{\rho}{t}
   +
   \der{\rho u_i}{x_i} \right) dV
   =
   0

using the Gauss theorem.

Since this relation should be satisfied everywhere (for any control volume), I request the integrand to be zero:

.. math::

   \der{\rho}{t}
   +
   \der{\rho u_i}{x_i}
   =
   0.

In this project, I assume that the liquid density is constant across the domain, i.e.

.. math::

   \der{\rho}{t}
   =
   0,

   \der{\rho}{x_i}
   =
   0_i.

Using these relations, I am left with the incompressibility constraint

.. math::

   \der{u_i}{x_i}
   =
   0.

