
.. _example_conservation_q:

.. include:: /references.txt

####################################
Conservation of Quadratic Quantities
####################################

In this section, we check the conservations of the volume-integrated quadratic quantities.
See :ref:`the governing equations <governing_equations>` for the definitions of these quantities.

*************
Configuration
*************

First, we simulate for :math:`50` time units to get a random flow.
:math:`Ra` is set to an extremely high value to mimic an inviscid condition.
This is followed by another run for :math:`10` time units, restarted from the previous simulation without the buoyancy force, so that the quadratic quantities should be conserved.

To see the effect of the time step size, I consider four different :ref:`Courant numbers <implicit_treatment>` :math:`0.1`, :math:`0.2`, :math:`0.4`, :math:`0.8`.

*******
Results
*******

As derived in :ref:`the governing equaations <governing_equations>`, the quadratic quantities should be constant over time.
This is not the case in this project, since we adopt an :ref:`explicit Runge-Kutta scheme <time_marchers>` to integrate the advective terms in time, which is dissipative (see |MORINISHI1998|, |COPPOLA2019|).
Thus these two quantities decrease monotonically over time, which is unfavourable but advantageous from a stability point of view.

The following graphs show the quadratic quantities as a function of time:

2D:

.. image:: https://raw.githubusercontent.com/NaokiHori/SimpleNSSolver/artifacts/artifacts/energy-2d/energy1.png
  :width: 600

3D:

.. image:: https://raw.githubusercontent.com/NaokiHori/SimpleNSSolver/artifacts/artifacts/energy-3d/energy1.png
  :width: 600

Here four different time step sizes are considered.
Also, the decays of :math:`K` and :math:`H` (:math:`t = 50` and :math:`60` are compared) are displayed.

2D:

.. image:: https://raw.githubusercontent.com/NaokiHori/SimpleNSSolver/artifacts/artifacts/energy-2d/energy2.png
  :width: 600

3D:

.. image:: https://raw.githubusercontent.com/NaokiHori/SimpleNSSolver/artifacts/artifacts/energy-3d/energy2.png
  :width: 600

Here the third-order convergence is observed, which is also reported by the previous works (e.g., |MORINISHI1998|, |HAM2002| and |COPPOLA2019|).
This indicates that the decay, which is a numerical artifact, is properly converged.

