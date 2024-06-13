
.. _example_nu_agreement:

#########################
Nusselt Number Agreements
#########################

*************
Configuration
*************

I consider three Prandtl numbers: :math:`Pr = 10^{-1}, 10^0, 10^1`.
The other parameters are listed below.

******
Result
******

:math:`Nu` calculated by four different formulae are considered:

* reference: heat flux on the walls

* red: energy injection

* blue: kinetic energy dissipation

* green: thermal energy dissipation

Their deviations from the reference value :math:`Nu_{wall}` are plotted to highlight the difference.
The reason why the lines are partly discontinuous is that the discrepancies are smaller than :math:`10^{-15}`.

2D, :math:`Pr = 10^{-1}`:

.. image:: https://raw.githubusercontent.com/NaokiHori/SimpleNSSolver/artifacts/artifacts/nusselt-2d-1.e-1/nusselt.png
  :width: 600

2D, :math:`Pr = 10^{ 0}`:

.. image:: https://raw.githubusercontent.com/NaokiHori/SimpleNSSolver/artifacts/artifacts/nusselt-2d-1.e+0/nusselt.png
  :width: 600

2D, :math:`Pr = 10^{ 1}`:

.. image:: https://raw.githubusercontent.com/NaokiHori/SimpleNSSolver/artifacts/artifacts/nusselt-2d-1.e+1/nusselt.png
  :width: 600

3D, :math:`Pr = 10^{-1}`:

.. image:: https://raw.githubusercontent.com/NaokiHori/SimpleNSSolver/artifacts/artifacts/nusselt-3d-1.e-1/nusselt.png
  :width: 600

3D, :math:`Pr = 10^{ 0}`:

.. image:: https://raw.githubusercontent.com/NaokiHori/SimpleNSSolver/artifacts/artifacts/nusselt-3d-1.e+0/nusselt.png
  :width: 600

3D, :math:`Pr = 10^{ 1}`:

.. image:: https://raw.githubusercontent.com/NaokiHori/SimpleNSSolver/artifacts/artifacts/nusselt-3d-1.e+1/nusselt.png
  :width: 600

The deviations should be small enough (around the rounding error).

