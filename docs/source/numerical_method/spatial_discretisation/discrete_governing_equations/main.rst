
.. _discrete_governing_equations:

.. include:: /references.txt

############################
Discrete governing equations
############################

************
Introduction
************

:ref:`The SMAC method <smac_method>` adopted in :ref:`the temporal discretisation <temporal_discretisation>` is a *de facto standard* to integrate the momentum balance in time while keeping the incompressibility constraint.
Spatial treatments are, on the other hand, more flexible and a *proper* scheme is totally case-dependent.
In this section, I focus on the following two features and derive the *discrete* governing equations :ref:`which are implemented in this project <implementation>`.

* Being simple

   In this project, I only consider the second-order-accurate central-difference schemes, where three points in one direction (one and its two neighbouring points) are used to evaluate the differentiations.
   Although its accuracy is not at all comparable to the spectral methods, it is versatile and robust.

* Mimicking the equations in the continuous domain

   Since I discretise and approximate the equations, it is impossible to keep all properties which the original equations have in the continuous domain.
   However, as reviewed in :ref:`the governing equations <governing_equations>`, there are several important relations among others, e.g.

   * The conservative terms do not alter the net amount of quantity, namely :math:`K \equiv \int k dV` and :math:`H \equiv \int h dV` are conserved in the absence of the dissipations.

   * The dissipative terms in the equations of the squared quantity (:math:`k` and :math:`h`) reduce the total amount of quantity.

   * Multiple Nusselt number descriptions yield statistically equivalent results.

   I focus on keeping these properties in the discrete domain.

.. note::

   In short, the conclusive treatment for the advective terms here is known as the energy-conserving scheme.
   To derive it, the equations are usually once mapped to a general coordinate system where the grid sizes in all directions are equidistant, which is followed by the discretisation (e.g. |KAJISHIMA1999|).

   Although it is a consistent and a powerful way, it tends to be fairly complicated (see the appendix, where I do that).
   In this section, I try to deduce the same conclusion without using this transform but by focusing on some simple relations of the discrete derivatives and integrals.
   Also the same concept is adopted to the internal energy equation.

***********
Derivations
***********

.. note::

   Hereafter I use notations like

   .. math::

      \sum_{\pic} q

   indicating that the quantity :math:`q`, which is defined at the cell-center, is summed in the :math:`x` direction.

   Similarly I define

   .. math::

      \sum_{\xic} q

   to imply that the quantity :math:`q`, which is defined at the cell-face, is summed in the :math:`x` direction.

.. toctree::
   :caption: Details and derivations
   :maxdepth: 1

   derivation/basic_operators/main.rst
   derivation/incompressibility_and_poisson_equation
   derivation/momentum_balance/main
   derivation/internal_energy_balance/main
   derivation/nusselt/main

.. toctree::
   :caption: Appendix
   :maxdepth: 1

   appendix/inconsistent_results/main
   appendix/strong_conservation_form/main

**********
Conclusion
**********

The following equations are implemented in the code, which are designed to keep the properties of the original governing equations discussed above.
See :ref:`the basic operators <basic_operators>` for details about the used symbols.

.. toctree::
   :maxdepth: 1

   conclusion/incompressibility
   conclusion/x_momentum/main
   conclusion/y_momentum/main
   conclusion/z_momentum/main
   conclusion/internal_energy/main
   conclusion/squared_velocity
   conclusion/squared_temperature

