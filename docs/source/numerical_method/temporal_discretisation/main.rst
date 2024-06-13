
.. _temporal_discretization:

.. include:: /references.txt

#######################
Temporal Discretization
#######################

In this section, we discuss how :ref:`the governing equations <governing_equations>` are integrated over time.

First, we briefly introduce several explicit and implicit schemes for integrating general ordinary differential equations over time.

We then elaborate on the implicit treatment applied to the diffusive terms to avoid severe restrictions on time-step sizes.

These fundamental techniques are first applied to the internal energy equation (evolution of temperature) due to its simplicity, and are then extended to the integration of the velocity field, where the momentum balance is integrated while keeping the incompressibility constraint.

.. toctree::
   :maxdepth: 1

   time_marcher
   implicit
   temperature
   momentum

