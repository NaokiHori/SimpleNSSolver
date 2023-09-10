
.. _statistics:

###################
`src/statistics.c`_
###################

.. _src/statistics.c: https://github.com/NaokiHori/SimpleNSSolver/blob/main/src/statistics.c

This file contains some functions which are used to collect temporally-averaged statistics and to output the results.

********
Overview
********

:math:`\ave{q}{t}`, which is a temporally-averaged value of a scalar field :math:`q`:

.. math::

   \ave{q}{t}
   =
   \frac{\int_{t} q dt}{\int_{t} dt},

is numerically approximated as:

.. math::

   \frac{\sum_{n=0}^{N-1} q_n}{N},

namely :math:`q` at :math:`N` different time steps are averaged.

By default, I collect the following quantities:

.. list-table:: Quantities
   :widths: 25 25 50
   :header-rows: 1

   * - Title
     - Position
     - Description
   * - :math:`\sum_{t} \ux`
     - :math:`x` cell face
     - N/A
   * - :math:`\sum_{t} \ux^2`
     - :math:`x` cell face
     - N/A
   * - :math:`\sum_{t} \uy`
     - :math:`y` cell face
     - N/A
   * - :math:`\sum_{t} \uy^2`
     - :math:`y` cell face
     - N/A
   * - :math:`\sum_{t} \uz`
     - :math:`z` cell face
     - 3D only
   * - :math:`\sum_{t} \uz^2`
     - :math:`z` cell face
     - 3D only
   * - :math:`\sum_{t} T`
     - cell center
     - N/A
   * - :math:`\sum_{t} T^2`
     - cell center
     - N/A
   * - :math:`\sum_{t} \ux T`
     - :math:`x` cell face
     - N/A

**************
Implementation
**************

A structure ``statistics_t``, which is defined in `include/statistics.h <https://github.com/NaokiHori/SimpleNSSolver/blob/main/include/statistics.h>`_, contains methods which collect and save the statistics.

.. note::

   Although each statistical data is collected as a two-dimensional (or three-dimensional) array, it is saved after summed in the homogeneous directions to reduce the storage requirement.
   This behaviour can be changed by a flag ``g_reduction`` defined in ``statistics.c``.

