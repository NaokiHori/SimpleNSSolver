
.. _fluid_compute_potential:

##############################
`src/fluid/compute_potential`_
##############################

.. _src/fluid/compute_potential: https://github.com/NaokiHori/SimpleNSSolver/tree/main/src/fluid/compute_potential

This directory contains functions to compute a scalar potential :math:`\psi` which projects a non-solenoidal velocity field to the solenoidal one by solving a Poisson equation

.. math::

   \frac{\delta^2 \psi}{\delta x_i \delta x_i}
   =
   \frac{1}{\gamma \Delta t} \frac{\delta u_i^*}{\delta x_i}.

The reason why I need to solve this equation is discussed in :ref:`the SMAC method <smac_method>`.

The mathematical background is briefly described in the following pages:

.. toctree::
   :maxdepth: 1

   math/governing_equation
   math/orthogonal_decomposition

.. mydeclare:: /../../src/fluid/compute_potential/dft.c
   :language: c
   :tag: fluid_compute_potential_dft

.. mydeclare:: /../../src/fluid/compute_potential/dct.c
   :language: c
   :tag: fluid_compute_potential_dct

.. warning::

   Notice the difference between the global array sizes ``glisize``, ``gljsize`` (given as the input values) and the local array sizes ``myisize``, ``myjsize`` (``pencil`` sizes).
   They are frequently switched in the description below because of the extensive use of the pencil rotations.

This function solves the Poisson equation described above.
Although it is a simple Poisson equation, the implementation is slightly complicated.
To make things easier, I encapsulate many complicated stuffs into one structure, which needs to be initialised when this function is called for the first time:

To solve the equation, I need the right-hand side of the equation, which is computed and stored in a buffer:

.. myliteralinclude:: /../../src/fluid/compute_potential/dft.c
   :language: c
   :tag: compute right-hand side of Poisson equation

.. myliteralinclude:: /../../src/fluid/compute_potential/dct.c
   :language: c
   :tag: compute right-hand side of Poisson equation

.. note::

   Two buffers are used, ``buf0`` and ``buf1``.
   The right-hand-side values are assigned to ``buf0``.

In this project, I use two algorithms to solve the Poisson equation.
One is a normal algorithm which depends on the discrete Fourier transforms (``DFT`` version), which is versatile and can be used anytime.

Another one is more efficient thanks to the discrete cosine transforms (``DCT`` version) and is approximately two times faster since the number of ``all-to-all`` communications are roughly halved.
However, this version is only applicable and invoked when the grid points in the :math:`x` direction are equidistantly placed.

.. toctree::

   dft
   dct

