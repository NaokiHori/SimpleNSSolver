
.. _logging_check_divergence:

##########
Divergence
##########

.. mydeclare:: /../../src/logging/divergence.c
   :language: c
   :tag: logging_check_divergence

This function computes the divergence of the velocity field and write the result to a file as well as the simulation time units at the moment.

Local divergence of the velocity field is computed as

.. math::

   \dder{\ux}{x}
   +
   \dder{\uy}{y}
   +
   \dder{\uz}{z}
   =
   \frac{
      \vat{\ux}{\pip, \pjc, \pkc}
      -
      \vat{\ux}{\pim, \pjc, \pkc}
   }{\Delta x_{\pic}}
   +
   \frac{
      \vat{\uy}{\pic, \pjp, \pkc}
      -
      \vat{\uy}{\pic, \pjm, \pkc}
   }{\Delta y}
   +
   \frac{
      \vat{\uz}{\pic, \pjc, \pkp}
      -
      \vat{\uz}{\pic, \pjc, \pkm}
   }{\Delta z},

which is defined at each cell center :math:`\left( \pic, \pjc, \pkc \right)`.

.. myliteralinclude:: /../../src/logging/divergence.c
   :language: c
   :tag: compute local divergence

Among others, the maximum value of the local divergence ``divmax``, which should be comparable to the machine epsilon (depending on the grid configuration), is checked:

.. myliteralinclude:: /../../src/logging/divergence.c
   :language: c
   :tag: check maximum

