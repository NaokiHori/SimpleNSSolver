
.. _spatial_discretization:

.. include:: /references.txt

######################
Spatial Discretization
######################

In this section, we extensively discuss spatial discretizations.

First, we describe how variables (coordinates, velocities, pressure, and temperature) are located within the computational domain using staggered grids.
We also briefly explain the domain decomposition technique used for parallelizing the domain.

Second, we address the numerical handling of equations, focusing on ambiguities in interpolations.
To achieve proper discretizations, we consider a computational coordinate system with uniform grid spacings, using equations in general rectilinear orthogonal coordinates.
We provide appropriate discretizations for the :ref:`governing equations <governing_equations>`: incompressibility, momentum balance, and internal energy balance.

Third, we prove that the resulting relations satisfy the desired properties in terms of the quadratic quantities :math:`k` and :math:`h`.
Specifically, the advective terms conserve the net :math:`k` and :math:`h`, while the diffusive terms conduct and dissipate them, as derived :ref:`here <quadratic_quantity_balance>`.
Additionally, we focus on the global energy balance.
For the squared velocity and the squared temperature, the source and sink terms must exactly balance in a statistically steady manner, which is not achieved unless energy consistency is fulfilled.
We show that with our proper treatment, all relations are indeed satisfied numerically.

Finally, the squared velocity and temperature are linked via the Nusselt number: a non-dimensional number measuring the heat transfer enhancement due to the convective effects.
There are several ways to calculate the Nusselt number, which all give an identical result in theory.
This relation is, not necessarily preserved from a numerical standpoint (c.f., |OSTILLAMONICO2015|).
We show that this exact relation can be replicated with our energy-consistent treatment.

.. toctree::
    :maxdepth: 1

    domain_setup/main
    discrete_governing_equations/main
    quadratic_quantities/main
    nusselt

