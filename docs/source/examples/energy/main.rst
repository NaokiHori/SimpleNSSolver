
.. _example_conservation_q:

.. include:: /references.txt

####################################
Conservation of quadratic quantities
####################################

In this section, I check the conservations of the volume-integrated quadratic quantities :math:`K` and :math:`H`.
See :ref:`the governing equations <governing_equations>` for the definitions of these quantities.

.. mydetails:: Signatures

   .. literalinclude:: data/ci_2d.txt
      :language: text

   .. literalinclude:: data/ci_3d.txt
      :language: text

*************
Configuration
*************

First, I simulate for :math:`50` time units to get a random velocity field and a random temperature field.
:math:`Ra` is set to an extremely high value to mimic the inviscid condition.
This is followed by another run for :math:`10` time units, restarted from the previous simulation without the buoyancy force, so that the quadratic quantities should be conserved.

To see the effect of the time step size, I consider four different safety factors :math:`0.1`, :math:`0.2`, :math:`0.4`, :math:`0.8`, which are multiplied to the maximum time step size computed in :ref:`src/decide_dt.c <decide_dt>`.

.. mydetails:: Configurations

   .. literalinclude:: data/exec_2d.sh
      :language: sh

   .. literalinclude:: data/exec_3d.sh
      :language: sh

*******
Results
*******

As derived in :ref:`the discrete momentum balance <momentum_balance>` and :ref:`the discrete internal energy balance <internal_energy_balance>`, they should be constant in time.
This is not the case in this project, since I adopt an :ref:`explicit Runge-Kutta scheme <temporal_discretisation>` to integrate the advective terms in time, which is dissipative (see |MORINISHI1998|, |COPPOLA2019|).
Thus these two quantities should decrease monotonically over time, which is unfavourable but advantageous from a stability point of view.

The following graphs show the quadratic quantities as a function of time:

* 2D

   .. image:: data/energy1_2d.png
      :width: 600

* 3D

   .. image:: data/energy1_3d.png
      :width: 600

Here four different time step sizes are considered.
Also, the decays of :math:`K` and :math:`H` (:math:`t = 50` and :math:`60` are compared) are displayed.

* 2D

   .. image:: data/energy2_2d.png
      :width: 600

* 3D

   .. image:: data/energy2_3d.png
      :width: 600

.. mydetails:: Script

   .. literalinclude:: data/process.py
      :language: python
      :linenos:

Here the third-order convergence is observed, which is also reported by the previous works (e.g. |MORINISHI1998|, |HAM2002| and |COPPOLA2019|).
This indicates that the decay, which is a numerical artifact, is converged relatively quickly.

.. seealso::

   The conservations of :math:`K` and :math:`H` depend strongly on the discretisation of the advective terms.
   Also the third-order temporal convergence is not expected if the spatial treatment is inconsistent.
   See :ref:`here <inconsistent_results>` for details.

