
.. _logging_nusselt:

###############
Nusselt numbers
###############

.. mydeclare:: /../../src/logging/nusselt/main.c
   :language: c
   :tag: logging_check_nusselt

This function computes the instantaneous Nusselt number (heat transfer enhancement by the convective effects) :math:`Nu \left( t \right)` based on the various definitions and write the results to a file.
I compute :math:`Nu` in several ways, which should return the same result in a statistical sense.

* :ref:`Heat flux on the walls <nu_heat_flux_discrete>`

* :ref:`Kinetic energy injection <nu_kinetic_energy_injection>`

* :ref:`Kinetic energy dissipation <nu_kinetic_energy_dissipation>`

* :ref:`Thermal energy dissipation <nu_thermal_energy_dissipation>`

