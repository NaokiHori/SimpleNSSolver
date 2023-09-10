
.. _eq_mass_discrete:

#############################
:ref:`Mass balance <eq_mass>`
#############################

Defined at :math:`\left( \pic, \pjc, \pkc \right)`:

.. math::

   \dder{\ux}{x}
   +
   \dder{\uy}{y}
   +
   \dder{\uz}{z}
   =
   0,

and thus

.. myliteralinclude:: /../../src/logging/divergence.c
   :language: c
   :tag: compute local divergence

.. seealso::

   * :ref:`fluid/correct_velocity <fluid_correct_velocity>`

   * :ref:`fluid/compute_potential <fluid_compute_potential>`

