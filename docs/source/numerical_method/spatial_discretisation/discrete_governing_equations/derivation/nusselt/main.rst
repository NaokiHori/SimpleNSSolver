
.. _nusselt_number_relations:

########################
Nusselt number relations
########################

In this part, I consider the discrete Nusselt number relations, whose continuous version is derived in :ref:`the governing equations <governing_equations>`.
In the continuous space, I first define the (temporally-averaged) heat flux, which is followed by the definition of the Nusselt number :math:`Nu`.
Then I derive other several relations to compute :math:`Nu`, all of which are equivalent to each other.
In this section, I apply the same procedure to derive the discrete analogues.

.. note::

   In general, the following relations assume that the system is in the statistically-steady state, i.e. they are satisfied after averaged in time :math:`\ave{q}{t}`.

.. toctree::
   :maxdepth: 1

   heat_flux
   kinetic_energy_injection
   kinetic_energy_dissipation
   thermal_energy_dissipation

.. seealso::

   These relations are used to monitor the simulations, which are implemented in :ref:`src/logging/nusselt <logging_nusselt>`.

