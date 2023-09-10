
.. _fluid_integrate:

######################
`src/fluid/integrate`_
######################

.. _src/fluid/integrate: https://github.com/NaokiHori/SimpleNSSolver/tree/main/src/fluid/integrate

This directory contains functions to predict the new scalar fields: the velocity and the temperature at the next Runge-Kutta level.
Although updating the temperature field is completed here, another step is needed for the velocity field since the updated field in general does not satisfy the divergence-free condition.

There are mainly three steps:

#. Compute right-hand-side terms of the Runge-Kutta scheme

   In general the update processes can be written as

   .. math::

      q^{k+1}
      =
      q^{k}
      +
      \Delta q.

   Here :math:`\Delta q` is computed and stored.

   .. note::

      At this point, only the advective, diffusive, and the pressure gradient effects are added to the buffer, where I have not taken into account the buoyancy force yet.

#. Couple external effects

   External effects such as the buoyancy force is appended here.
   The reason why I separate this process from the :ref:`fluid_compute_rhs <fluid_compute_rhs>` is to make the other couplings (e.g. driving pressure, no-slip effects induced by the immersed boundary method) easier by exposing the API.

#. Predict field

   :math:`q` is updated from :math:`q^k` to :math:`q^{k+1}` by solving an ordinary differential equation for each quantity.

   .. note::

      While the temperature field is fully updated by this function, the momentum field is not divergence-free at this point, which will be considered in :ref:`the next step <fluid_compute_potential>`.

.. toctree::
   :maxdepth: 1

   compute_rhs
   predict_field

