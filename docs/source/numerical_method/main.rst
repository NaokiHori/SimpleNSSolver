
.. _numerics:

################
Numerical method
################

Numerical treatments of :ref:`the governing equations <governing_equations>` are discussed here.

* :ref:`Temporal discretisation <temporal_discretisation>`

   Conventional techniques to integrate the ordinary differential equations in time (e.g. SMAC method, implicit treatment of the diffusive terms) are briefly discussed.

* :ref:`Spatial discretisation <spatial_discretisation>`

   How the variables are discretised in the computational (discrete) domain is discussed.
   Although there are :ref:`five equations <governing_equations>`, :math:`k` and :math:`h` are not solved numerically because they depend on the other three.

   These *implicit* dependencies are easily broken if I do not treat the three main equations properly.
   For instance, the kinetic energy (squared velocity) is artificially added or removed if the momentum equation is not properly discretised in space.
   This correct numerical treatment is of my main interest in this section.

.. toctree::
   :hidden:

   temporal_discretisation/main
   spatial_discretisation/main

