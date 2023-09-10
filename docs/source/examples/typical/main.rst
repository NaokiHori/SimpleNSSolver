
.. _example_typical_case:

.. include:: /references.txt

#############################
Typical 2D and 3D simulations
#############################

.. mydetails:: Signatures

   .. literalinclude:: data/ci_2d.txt
      :language: text

   .. literalinclude:: data/ci_3d.txt
      :language: text

*************
Configuration
*************

.. mydetails:: Two-dimensional case:

   .. literalinclude:: data/exec_2d.sh
      :language: sh
      :linenos:

.. mydetails:: Three-dimensional case:

   .. literalinclude:: data/exec_3d.sh
      :language: sh
      :linenos:

**************
Visualisations
**************

Temperature fields at the end of the simulations are visualised here.

* 2D:

   .. image:: data/snapshot_2d.png
      :width: 600

   .. mydetails:: Script

      .. literalinclude:: data/snapshot_2d.py
         :language: python
         :linenos:

* 3D

   .. image:: data/snapshot_3d.png
      :width: 600

   .. mydetails:: Script

      .. literalinclude:: data/snapshot_3d.py
         :language: python
         :linenos:

****************************
Incompressibility constraint
****************************

The local divergence of the velocity fields should be sufficiently small, which is monitored during the simulation and checked here.

* 2D

   .. image:: data/divergence_2d.png
      :width: 600

* 3D

   .. image:: data/divergence_3d.png
      :width: 600

.. mydetails:: Script

   .. literalinclude:: data/divergence.py
      :language: python
      :linenos:

.. seealso::

   :ref:`src/logging/divergence.c <logging_check_divergence>`

***************
Nusselt numbers
***************

=====================
As a function of time
=====================

:math:`Nu` calculated using the different formulae, which are monitored during the run, are shown as a function of time:

   * red: heat fluxes on the walls

   * blue: energy input

   * green: kinetic energy dissipation

   * magenta: thermal energy dissipation

* 2D

   .. image:: data/nusselt_time_2d.png
      :width: 600

* 3D

   .. image:: data/nusselt_time_3d.png
      :width: 600

.. mydetails:: Script

   .. literalinclude:: data/nusselt_time.py
      :language: python
      :linenos:

.. note::

   The black-dashed line in the two-dimensional result shows a reference value by |VANDERPOEL2013| with the same :math:`Ra` and :math:`Pr` but the different domain geometry is different (box).

.. seealso::

   Check :ref:`Nusselt number relations <nusselt_number_relations>` to see the definition of each :math:`Nu`, and :ref:`src/logging/nusselt <logging_nusselt>` to see how they are computed.

=======
Average
=======

There are two contributions, advective contribution:

.. math::

   \ave{\ux T}{y,z,t},

and diffusive contribution:

.. math::

   - \frac{1}{\sqrt{Pr} \sqrt{Ra}} \frac{d}{dx} \ave{T}{y,z,t},

which are a function of the wall-normal position :math:`x` and shown separately:

* 2D

   .. image:: data/nusselt_x_2d.png
      :width: 600

* 3D

   .. image:: data/nusselt_x_3d.png
      :width: 600

.. mydetails:: Script

   .. literalinclude:: data/nusselt_x.py
      :language: python
      :linenos:

*******************
Standard deviations
*******************

Variances of (red) :math:`\ux`, (blue) :math:`\uy`, (magenta) :math:`\uz`, and (green) :math:`T` are shown here.

* 2D

   .. image:: data/std_2d.png
      :width: 600

* 3D

   .. image:: data/std_3d.png
      :width: 600

.. mydetails:: Script

   .. literalinclude:: data/std.py
      :language: python
      :linenos:

.. note::

   Although the :math:`y` and the :math:`z` directions are homogeneous, the blue and the magenta lines may show a deviation.
   This is attributed to the low :math:`Ra` and the short simulation time to reduce the running cost.

