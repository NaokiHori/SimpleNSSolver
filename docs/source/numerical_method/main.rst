
.. _numerics:

################
Numerical Method
################

Numerical treatments of :ref:`the governing equations <governing_equations>` are discussed here.

In the first part, the discretization of the variables in the computational (discrete) domain is discussed.
Although there are :ref:`five equations <governing_equations>`, the relations with respect to :math:`k` and :math:`h` are not solved numerically because they depend on the other three.
These dependencies, which are obvious in theory, are so fragile that they are not satisfied if we do not treat the main three relations properly.
The proper spatial numerical treatment is our main interest.

In the second part, conventional techniques to integrate the given ordinary differential equations in time are briefly discussed.
Specifically, we focus on the SMAC method for dealing with the momentum balance and the incompressibility constraint, and the implicit treatment of the diffusive terms.

To update the flow field, we need to solve linear systems (when diffusive terms are treated implicitly) and a Poisson equation, which are separately discussed.
Both rely on solving tri-diagonal matrices, which is also elaborated.

.. toctree::
   :maxdepth: 1

   spatial_discretisation/main
   temporal_discretisation/main
   linear_system
   poisson
   tdm

