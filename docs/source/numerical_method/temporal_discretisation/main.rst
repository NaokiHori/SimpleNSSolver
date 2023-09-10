
.. include:: /references.txt

.. _temporal_discretisation:

#######################
Temporal discretisation
#######################

In this section, I discuss how :ref:`the governing equations <governing_equations>` are integrated in time.

To begin with, general remarks (the overall time marching schemes and the implicit treatment of the diffusive terms) are discussed.
The treatment of the temperature field is presented in the first part as it is simpler.
Special attention is needed for the conservation of the mass and the balance of the momentum, which is discussed in the second part.

.. toctree::
   :maxdepth: 1

   time_marcher
   implicit
   temperature
   smac_method

