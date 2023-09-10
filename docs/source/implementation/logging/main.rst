
.. _logging:

##############
`src/logging`_
##############

.. _src/logging: https://github.com/NaokiHori/SimpleNSSolver/tree/main/src/logging

This directory contains source files to compute and output quantities to be monitored during the simulation.

.. mydeclare:: /../../src/logging/main.c
   :language: c
   :tag: check_and_output

This function is a main function and an entry point of this directory, which calls other functions which compute and output log data discussed below.

.. mydeclare:: /../../src/logging/main.c
   :language: c
   :tag: show_progress

This function displays the current status such as

   * time step ``step``

   * simulation time ``time``

   * time step size ``dt``

   * elapsed time ``wtime``

.. myliteralinclude:: /../../src/logging/main.c
  :language: c
  :tag: show progress to standard output and file

Note that this information is displayed to ``stdout`` and is also written to a file (by default ``output/log/progress.dat``).

Monitored quantities and their implementations are described below:

.. toctree::
   :maxdepth: 1

   divergence
   momentum
   energy
   nusselt

