
.. _example_typical_case:

.. include:: /references.txt

#############
Typical Cases
#############

**************
Visualisations
**************

Temperature fields at the end of the simulations.

2D:

.. image:: https://raw.githubusercontent.com/NaokiHori/SimpleNSSolver/artifacts/artifacts/typical-2d/snapshot.png
  :width: 600

3D:

.. image:: https://raw.githubusercontent.com/NaokiHori/SimpleNSSolver/artifacts/artifacts/typical-3d/snapshot.png
  :width: 600

****************************
Incompressibility constraint
****************************

Maximum divergence of the velocity field, which should be sufficiently small.

2D:

.. image:: https://raw.githubusercontent.com/NaokiHori/SimpleNSSolver/artifacts/artifacts/typical-2d/divergence.png
  :width: 600

3D:

.. image:: https://raw.githubusercontent.com/NaokiHori/SimpleNSSolver/artifacts/artifacts/typical-3d/divergence.png
  :width: 600

***************
Nusselt numbers
***************

=========
Evolution
=========

:math:`Nu` calculated using the different formulae, which are monitored during the run, are shown as a function of time:

* red: heat fluxes on the walls

* blue: energy input

* green: kinetic energy dissipation

* magenta: thermal energy dissipation

2D:

.. image:: https://raw.githubusercontent.com/NaokiHori/SimpleNSSolver/artifacts/artifacts/typical-2d/nusselt_time.png
  :width: 600

3D:

.. image:: https://raw.githubusercontent.com/NaokiHori/SimpleNSSolver/artifacts/artifacts/typical-3d/nusselt_time.png
  :width: 600

.. note::

   The black-dashed line in the two-dimensional result shows a reference value by |VANDERPOEL2013| with the same :math:`Ra` and :math:`Pr` but the different domain geometry is different (box).

.. seealso::

   :ref:`Nusselt number relations <discrete_nusselt>`.

=========================
Temporary-Averaged Values
=========================

As derived :ref:`here <eq_heat_transfer>`, there are two contributions which transfer heat: advective contribution:

.. math::

    \sumzc
    \sumyc
    \frac{J}{\sfact{1}}
    \vel{1}
    \ave{T}{\gcs{1}},

and diffusive contribution:

.. math::

    -
    \sumzc
    \sumyc
    \frac{1}{\sqrt{Pr} \sqrt{Ra}}
    \frac{J}{\sfact{1}}
    \frac{1}{\sfact{1}}
    \dif{T}{\gcs{1}}.

After averaged over time and homogeneous directions, they are displayed as a function of the wall-normal position :math:`x` here:

2D:

.. image:: https://raw.githubusercontent.com/NaokiHori/SimpleNSSolver/artifacts/artifacts/typical-2d/nusselt_x.png
  :width: 600

3D:

.. image:: https://raw.githubusercontent.com/NaokiHori/SimpleNSSolver/artifacts/artifacts/typical-3d/nusselt_x.png
  :width: 600

*******************
Standard deviations
*******************

Variances of (red) :math:`\ux`, (blue) :math:`\uy`, (magenta) :math:`\uz`, and (green) :math:`T` are shown here.

2D:

.. image:: https://raw.githubusercontent.com/NaokiHori/SimpleNSSolver/artifacts/artifacts/typical-2d/std.png
  :width: 600

3D:

.. image:: https://raw.githubusercontent.com/NaokiHori/SimpleNSSolver/artifacts/artifacts/typical-3d/std.png
  :width: 600

.. note::

   Although the :math:`y` and the :math:`z` directions are homogeneous, the blue and magenta lines may deviate, which is attributed to the low :math:`Ra` and short time.

