
.. _fluid_init:

###################
`src/fluid/init.c`_
###################

.. _src/fluid/init.c: https://github.com/NaokiHori/SimpleNSSolver/blob/main/src/fluid/init.c

This section discusses how the velocity and the temperature fields are initialised.

.. mydeclare:: /../../src/fluid/init.c
   :language: c
   :tag: fluid_init

First, arrays associated to ``fluid_t`` (arrays to store the flow field, buffers to keep the intermediate information, etc.) are allocated:

.. myliteralinclude:: /../../src/fluid/init.c
   :language: c
   :tag: allocate arrays

Then the fields are loaded from the corresponding files which are contained in a specified directory, which is given as the argument (see ``exec.sh``):

.. myliteralinclude:: /../../src/fluid/init.c
   :language: c
   :tag: load flow fields

The boundary conditions are imposed and the halo cells are communicated:

.. myliteralinclude:: /../../src/fluid/init.c
   :language: c
   :tag: impose boundary conditions and communicate halo cells

.. note::

   The shape of the arrays (number of boundary cells and halo cells) is identical except for :math:`\ux`.

Also the momentum and the thermal diffusivities are computed and stored:

.. myliteralinclude:: /../../src/fluid/init.c
   :language: c
   :tag: compute diffusivities

