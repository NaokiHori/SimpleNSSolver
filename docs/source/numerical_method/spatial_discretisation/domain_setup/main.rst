
.. _domain_setup:

############
Domain setup
############

.. note::

   For simplicity I neglect the :math:`z` direction in this page.
   Note that the same condition applies as for the :math:`y` direction.

***************************************
Staggered grid and domain decomposition
***************************************

In this project, I use the staggered grid arrangement, i.e. the pressure and all velocity components are defined at different locations:

.. image:: images/staggered1.png
   :width: 800

where the left figure shows the locations and indices of the variables whose notations are used in the equations, while the descriptions found in the right figure are used in the code.

In the current project, I assume that the domain is wall-bounded in the :math:`x` direction, while periodic boundary conditions are imposed in the :math:`y` direction.
The spatial resolution is defined in terms of the number of cell centers (where the pressure and the temperature are defined): ``glisize = domain->glsizes[0]`` and ``gljsize = domain->glsizes[1]`` in the :math:`x` and the :math:`y` directions, respectively.
Here the prefix **gl** is used to emphasise they are **global** sizes.

This library supports the process parallelisation (in particular by means of ``MPI``), which is achieved by splitting the whole domain into smaller blocks, and each process is responsible for ``myisize = domain->mysizes[0]`` times ``myjsize = domain->mysizes[1]`` cell centers:

.. image:: images/domain.png
   :width: 800

.. note::

   * No decomposition in the :math:`x` direction

      In this project, the wall-normal direction is not decomposed (so-called `pencil decomposition <https://github.com/NaokiHori/SimpleDecomp>`_ is adopted).
      Thus ``glisize = glsizes[0]`` and ``myisize = mysizes[0]`` are equal.

   * Pencil sizes can be different

      Although the domain is decomposed so that each process has a similar workload, ``myjsize`` can be different for each process, epsecially when the domain size is not divisible by the number of processes in the direction.

In each process, each variable is positioned as follows:

* ``UX(i, j)``

   .. image:: images/staggered2.png
      :width: 600

* ``UY(i, j)``

   .. image:: images/staggered3.png
      :width: 600

* ``P(i, j)`` and ``T(i, j)``

   .. image:: images/staggered4.png
      :width: 600

Sizes of the two-dimensional arrays are summarised as follows:

============ =============== ===============
Name         ``i`` range     ``j`` range
============ =============== ===============
``UX``       ``1:myisize+1`` ``0:myjsize+1``
``UY``       ``0:myisize+1`` ``0:myjsize+1``
``P``, ``T`` ``0:myisize+1`` ``0:myjsize+1``
============ =============== ===============

.. note::

   * Cell-face positions and cell-center positions

      Although two different positions are adopted, only the cell faces are the free parameter since the cell centers are assumed to be positioned in the middle of the neighbouring cell faces.

   * Halo cells

      In order to evaluate the differentiations in the :math:`y` direction close to the domain edges, additional cells (halo cells) are attached in the :math:`y` boundaries.

   * Boundary points

      ``UY(0, j)``, ``UY(glisize+1, j)`` are shifted towards the near wall locations so that I can impose the velocity boundary conditions directly.
      Same holds for the temperature and the pressure.

      Another way is to locate a cell inside one wall, whose values are extrapolated from the interior values (ghost cells).
      It is easy to confirm that this extrapolation and the current implementation are identical.

*****************************************
Uniform and stretched grid configurations
*****************************************

In the :math:`y` direction, distance between two neighbouring points should be equal.
In the :math:`x` direction, non-uniform grid arrangement can be adopted, which could be useful to resolve boundary layers close to the walls.
To identify the positions, it is necessary to use two variables which are in ``domain`` structure, whose definitions are described below.

#. Cell-face positions ``domain->xf``

   ``XF(i)`` is used to describe the position of the cell faces in :math:`x` direction (**f** comes from face), i.e. the locations where :math:`\ux` is defined.

   .. note::

      This should be given as the initial condition.

   .. image:: images/grid1.png
      :width: 800

#. Cell-center positions ``domain->xc``

   ``XC(i)`` is used to describe the position of the cell centers in :math:`x` direction (**c** comes from center), i.e. the locations where :math:`p`, :math:`T`, and :math:`\uy` are defined.

   Note again that the first and the last points are defined at cell-faces.
   Except these two points, the following relation gives the relation between the cell faces and the center:

   .. math::

      XC \left( i \right)
      =
      \frac{1}{2} XF \left( i     \right)
      +
      \frac{1}{2} XF \left( i + 1 \right).

   .. image:: images/grid2.png
      :width: 800

   .. note::

      Although this is not a free parameter, this should also be given as the initial condition as well as ``xf``.

#. Cell-face distance ``domain->dxf``

   ``DXF(i)`` is used to describe the distance between two neighbouring cell faces:

   .. math::

      DXF \left( i \right)
      =
      XF \left( i + 1 \right)
      -
      XF \left( i     \right).

   .. image:: images/grid3.png
      :width: 800

   .. note::

      This is computed in :ref:`src/domain/init.c <domain>`.

#. Cell-center distance ``domain->dxc``

   ``DXC(i)`` is used to describe the distance between two neighbouring cell centers:

   .. math::

      DXC \left( i \right)
      =
      XC \left( i     \right)
      -
      XC \left( i - 1 \right).

   .. image:: images/grid4.png
      :width: 800

   .. note::

      This is computed in :ref:`src/domain/init.c <domain>`.

