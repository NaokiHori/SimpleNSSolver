
.. _example_nu_agreement:

#########################
Nusselt number agreements
#########################

Here, I investigate the perfect :math:`Nu` matches.
To eliminate oscillations, relatively low Rayleigh numbers are considered.
In particular a lower Rayleigh number is employed for the three-dimensional cases to stabilise the thermal plumes.

.. mydetails:: Signatures

   .. literalinclude:: data/ci_1.e-1_2d.txt
      :language: text

   .. literalinclude:: data/ci_1.e+0_2d.txt
      :language: text

   .. literalinclude:: data/ci_1.e+1_2d.txt
      :language: text

   .. literalinclude:: data/ci_1.e-1_3d.txt
      :language: text

   .. literalinclude:: data/ci_1.e+0_3d.txt
      :language: text

   .. literalinclude:: data/ci_1.e+1_3d.txt
      :language: text

*************
Configuration
*************

I consider three Prandtl numbers: :math:`Pr = 10^{-1}, 10^0, 10^1`.
The other parameters are listed below.

.. mydetails:: Details

   .. literalinclude:: data/exec_2d.sh
      :language: sh

   .. literalinclude:: data/exec_3d.sh
      :language: sh

******
Result
******

:math:`Nu` calculated by four different formulae are considered:

   * reference: heat flux on the walls

   * red: energy injection

   * blue: kinetic energy dissipation

   * green: thermal energy dissipation

.. seealso::

   :ref:`Nusselt number relations <nusselt_number_relations>` for the derivations.

Their deviations from the reference value :math:`Nu_{wall}` are plotted to highlight the difference.
The reason why the lines are partly discontinuous is that the discrepancies are smaller than :math:`10^{-15}`.

* 2D, :math:`Pr = 10^{-1}`

   .. image:: data/nusselt_1.e-1_2d.png
      :width: 600

* 2D, :math:`Pr = 10^{ 0}`

   .. image:: data/nusselt_1.e+0_2d.png
      :width: 600

* 2D, :math:`Pr = 10^{ 1}`

   .. image:: data/nusselt_1.e+1_2d.png
      :width: 600

* 3D, :math:`Pr = 10^{-1}`

   .. image:: data/nusselt_1.e-1_3d.png
      :width: 600

* 3D, :math:`Pr = 10^{ 0}`

   .. image:: data/nusselt_1.e+0_3d.png
      :width: 600

* 3D, :math:`Pr = 10^{ 1}`

   .. image:: data/nusselt_1.e+1_3d.png
      :width: 600

The deviations should be small enough (around the rounding error).

.. mydetails:: Script

   .. literalinclude:: data/process.py
      :language: python
      :linenos:

.. seealso::

   The :math:`Nu` consistency discussed here depends to a large extent on how I compute the dissipations.
   See a more detailed analysis :ref:`here <inconsistent_results>`.

