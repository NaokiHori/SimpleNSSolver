
.. _discrete_governing_equations:

.. include:: /references.txt

############################
Discrete Governing Equations
############################

:ref:`The SMAC method <momentum_integration>` adopted in :ref:`the temporal discretization <temporal_discretization>` is a *de facto standard* to integrate the momentum balance in time while keeping the incompressibility constraint.
Spatial treatments are, on the other hand, more flexible, and a proper scheme is totally case-dependent.
In this section, we focus on the following two features and derive the discrete governing equations which are implemented in this project.

* Being simple

    In this project, we only consider second-order-accurate central-difference schemes, where three points in one direction (one and its two neighboring points) are used to evaluate the differentiations.
    Although its accuracy is not at all comparable to spectral methods, it is versatile, robust, and also applicable to multiphase flows due to its locality.

* Preserving properties of the original equations

    Since we approximate the equations by discretizations, it is impossible to keep all properties that the original equations have in the continuous domain.
    However, as reviewed in :ref:`the governing equations <governing_equations>`, there are several important relations, among others, e.g.,

    * The conservative terms do not alter the net amount of quadratic quantities, namely :math:`\int k dV` and :math:`\int h dV` are conserved in the absence of viscous effects.

    * The dissipative terms in the equations of the squared quantities (:math:`k` and :math:`h`) reduce the total amount of the quantity.

    We focus on replicating these properties in the discrete domain.

In short, the conclusive treatment for the advective terms here is known as the energy-conserving scheme, and we refer to the overall scheme which treat all other terms in a consistent manner as *energy-consistent* scheme.

.. toctree::
    :maxdepth: 1

    strong_conservation_form
    resulting_schemes/main

.. seealso::

    * |KAJISHIMA1999|
    * |KAJISHIMA2017|

