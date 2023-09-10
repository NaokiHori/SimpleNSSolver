
.. _poisson_implementation:

##############
Implementation
##############

In :ref:`the previous section <poisson_eigen_decomposition>`, I derive a discrete Poisson equation in a simplified form using the eigen decomposition in the :math:`y` direction.
Here, I consider the resulting systems in detail and see the corresponding implementations.

By letting

.. math::

   \vat{P}{\pic, k_y}
   \equiv
   \mathcal{F}_d \left[ \vat{p}{\pic, \pjc} \right],

I am now left with

.. math::

   \frac{
      \frac{
         \vat{P}{\pipp, k_y}
         -
         \vat{P}{\pic,  k_y}
      }{
         \vat{\Delta x}{\pip}
      }
      -
      \frac{
         \vat{P}{\pic,  k_y}
         -
         \vat{P}{\pimm, k_y}
      }{
         \vat{\Delta x}{\pim}
      }
   }{
      \vat{\Delta x}{\pic}
   }
   -
   \frac{4}{\Delta y^2}
   \sin^2 \left( \frac{\pi}{N_y} k_y \right)
   \vat{P}{\pic, k_y}
   =
   \vat{Q}{\pic, k_y},

or equivalently

.. math::

   \vat{l}{\pic} \vat{P}{\pimm, k_y}
   +
   \vat{c}{\pic} \vat{P}{\pic,  k_y}
   +
   \vat{u}{\pic} \vat{P}{\pipp, k_y}
   =
   \vat{Q}{\pic, k_y},

where

.. math::

   \vat{l}{\pic}
   \equiv
   \frac{1}{\vat{\Delta x}{\pim} \vat{\Delta x}{\pic}},

.. math::

   \vat{u}{\pic}
   \equiv
   \frac{1}{\vat{\Delta x}{\pip} \vat{\Delta x}{\pic}},

and

.. math::

   \vat{c}{\pic}
   \equiv
   -
   \vat{l}{\pic}
   -
   \vat{u}{\pic}
   -
   \frac{4}{\Delta y^2}
   \sin^2 \left( \frac{\pi}{N_y} k_y \right),

which is a tri-diagonal linear system and can be solved in :math:`\mathcal{O} \left( N_x \right)` operations.

The lower- and upper-diagonal parts are assigned here:

.. myliteralinclude:: /../../src/fluid/compute_potential/dft.c
   :language: c
   :tag: initialise tri-diagonal matrix in x direction

.. myliteralinclude:: /../../src/fluid/compute_potential/dct.c
   :language: c
   :tag: initialise tri-diagonal matrix in y direction

.. myliteralinclude:: /../../src/fluid/compute_potential/dct.c
   :language: c
   :tag: initialise tri-diagonal matrix in z direction

which are initialised when the solver is called for the first time and are reused.

The center-diagonal part is computed every time when the system is solved:

.. myliteralinclude:: /../../src/fluid/compute_potential/dft.c
   :language: c
   :tag: set center diagonal components

.. myliteralinclude:: /../../src/fluid/compute_potential/dct.c
   :language: c
   :tag: set center diagonal components

The wave numbers are computed beforehand:

.. myliteralinclude:: /../../src/fluid/compute_potential/dft.c
   :language: c
   :tag: initialise eigenvalues in homogeneous directions

.. myliteralinclude:: /../../src/fluid/compute_potential/dct.c
   :language: c
   :tag: initialise eigenvalues in homogeneous directions

