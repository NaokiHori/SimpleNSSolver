
.. _domain:

###############
`src/domain.c`_
###############

.. _src/domain.c: https://github.com/NaokiHori/SimpleNSSolver/blob/main/src/domain.c

This file contains functions for manipulating a structure ``domain_t``, which stores spatial information about the domain as well as information about `the domain decomposition <https://github.com/NaokiHori/SimpleDecomp>`_.

**************
Initialisation
**************

.. mydeclare:: /../../src/domain.c
   :language: c
   :tag: domain_init

A structure ``domain_t`` is defined in `include/domain.h <https://github.com/NaokiHori/SimpleNSSolver/blob/main/include/domain.h>`_, which contains bunch of parameters:

.. myliteralinclude:: /../../include/domain.h
   :language: c
   :tag: definition of a structure domain_t

This function ``domain_init`` allocates and initialises a structure ``domain_t`` and its members, whose pointer is returned.

There are several important steps, which are summarised as follows:

#. Loading spatial information

   The size of the **global** array ``glsizes`` (see :ref:`the domain setup <domain_setup>`), which should be provided by the user as a set of initial conditions, is loaded.

   .. myliteralinclude:: /../../src/domain.c
      :language: c
      :tag: load spatial information

   At the same time, the physical sizes (lengths) of the domain (``lengths``) and the grid positions in the :math:`x` direction (``xf`` and ``xc``) are loaded.

   .. seealso::

      ``domain_save`` described below.

#. Compute grid sizes

   Based on the loaded information in the previous step, grid sizes in all directions are computed here:

   .. myliteralinclude:: /../../src/domain.c
      :language: c
      :tag: compute grid sizes

   Here ``dxf`` and ``dxc`` define face-to-face and center-to-center distances in the :math:`x` direction, respectively.
   In the other directions, since grids are equidistantly positioned, grid sizes are simply given by the length of the domain divided by the number of grids.

#. Decomposing the domain into ``pencils``.

   Since I know the global size of the system now, my next objective is to split the whole domain and to distribute each segment to each process, which is achieved by a function ``sdecomp.construct`` provided by `SimpleDecomp <https://github.com/NaokiHori/SimpleDecomp>`_ library:

   .. myliteralinclude:: /../../src/domain.c
      :language: c
      :tag: initialise sdecomp to distribute the domain

   .. note::

      In the three-dimensional case, the number of processes in the :math:`y` and the :math:`z` directions is adjusted to minimise the cost of pencil rotations, which are crucial for :ref:`solving the Poisson equations efficiently <fluid_compute_potential>`.

#. Compute local array sizes and offsets

   One the domain decomposition is complete, I can calculate the size of the **local** array ``mysizes`` and offsets in terms of the global positions ``offsets``:

   .. myliteralinclude:: /../../src/domain.c
      :language: c
      :tag: local array sizes and offsets

*************
Load and save
*************

.. mydeclare:: /../../src/domain.c
   :language: c
   :tag: domain_load

.. mydeclare:: /../../src/domain.c
   :language: c
   :tag: domain_save

These functions load or save some of the members in ``domain_t`` related to the domain sizes and the resolutions.

