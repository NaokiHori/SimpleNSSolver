
.. _halo:

#############
`src/halo.c`_
#############

.. _src/halo.c: https://github.com/NaokiHori/SimpleNSSolver/blob/main/src/halo.c

This file contains functions which take care of the inter-process communication in order to exchange the edge (halo) information.

********
Overview
********

#. Identifying the neighbours

   I need to check from which process (to which process) I receive (send) the halo information.
   This is achieved by using ``sdecomp.get_neighbours``.

#. Packing and unpacking

   Data to be communicated is not usually contiguous in memory.
   To encapsulate the packing-unpacking procedures, I define ``MPI_Datatype`` for each scalar field.
   Because the staggered grid is used, different data types should be used for each variable.

*******
Details
*******

.. mydetails:: :math:`y` halo cells

   .. include:: y-halo.rst

.. mydetails:: :math:`z` halo cells (3D only)

   .. include:: z-halo.rst

.. note::

   In spite of the redundancy and the increase in the size of the data, the boundary cells in the :math:`x` direction are also communicated for simplicity.

