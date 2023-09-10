########
Examples
########

This section presents several examples and their analyses.

To check the code quality regularly and automatically, this section is partially updated by `the GitHub Actions <https://github.com/features/actions>`_ when the ``main`` branch is updated.
The workflow file is `here <https://github.com/NaokiHori/SimpleNSSolver/blob/main/.github/workflows/ci.yml>`_.

.. note::

   To reduce the computational effort, I restrict the duration of the simulation ``timemax`` to relatively small values.
   As a result,

      * I might start to collect statistics even though the system has not reached the statistically-stationary state yet, and / or

      * statistics which are described in the following pages might not be converged.

.. toctree::
   :maxdepth: 1

   typical/main
   energy/main
   nu/main
   gl/main

