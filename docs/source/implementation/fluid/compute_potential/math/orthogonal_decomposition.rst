
.. _poisson_orthogonal_decomposition:

########################
Orthogonal decomposition
########################

By adopting the second-order-accurate central-difference approximation, the two-dimensional Poisson equation formulated in :ref:`the previous section <poisson_governing_equation>` is discretised as

.. math::

   \frac{
      \frac{
         \vat{p}{\pipp, \pjc}
         -
         \vat{p}{\pic,  \pjc}
      }{
         \vat{\Delta x}{\pip}
      }
      -
      \frac{
         \vat{p}{\pic,  \pjc}
         -
         \vat{p}{\pimm, \pjc}
      }{
         \vat{\Delta x}{\pim}
      }
   }{
      \vat{\Delta x}{\pic}
   }
   +
   \frac{
      \frac{
         \vat{p}{\pic, \pjpp}
         -
         \vat{p}{\pic, \pjc }
      }{
         \Delta y
      }
      -
      \frac{
         \vat{p}{\pic, \pjc }
         -
         \vat{p}{\pic, \pjmm}
      }{
         \Delta y
      }
   }{
      \Delta y
   }
   =
   \vat{q}{\pic, \pjc},

which is defined at each cell center.

Note that I assume the grids in the :math:`y` direction are uniformly-distributed, while those in the :math:`x` direction can vary.

Solving the above equation is equivalent to handle a linear system whose size is :math:`N_x \times N_y`, where :math:`N_x` and :math:`N_y` are the number of cell centers in each direction.
Solving it in a direct way such that the incompressibility is perfectly satisfied (in a computational sense) is impractical.

To solve the system efficiently, I consider the orthogonal decomposition of the system in the uniform direction (:math:`y`), i.e. I aim to reduce the scattered nature of the original linear system by diagonalising it in the :math:`y` direction.

For notational simplicity, I neglect the :math:`x` indices for a while.

Since the grid points are equidistantly spaced and a periodic boundary condition is assumed in the :math:`y` direction, I expand the equation using the trigonometric functions (Fourier series):

.. math::

   \vat{p}{\pjc}
   \equiv
   \sum_{l = 0}^{N_y - 1}
   \vat{P}{l}
   \exp \left( \frac{2 \pi}{N_y} I \pjc l \right),

where :math:`I` is the imaginary unit :math:`\sqrt{-1}`, while :math:`P_l` is the :math:`l`-th (discrete) wave number.
Note that :math:`P_l` satisfies the Hermitian symmetry since :math:`p_j \in \mathbb{R}`.

The two neighbouring points which are necessary to evaluate the discrete Poisson equation are given by

.. math::

   \vat{p}{\pjmm}
   &
   =
   \sum_{l = 0}^{N_y - 1}
   \vat{P}{l}
   \exp \left[ \frac{2 \pi}{N_y} I \left( \pjmm \right) l \right]

   &
   =
   \sum_{l = 0}^{N_y - 1}
   \vat{P}{l}
   \exp \left( \frac{2 \pi}{N_y} I \pjc l \right)
   \exp \left( - \frac{2 \pi}{N_y} I l \right),

   \vat{p}{\pjpp}
   &
   =
   \sum_{l = 0}^{N_y - 1}
   \vat{P}{l}
   \exp \left[ \frac{2 \pi}{N_y} I \left( \pjpp \right) l \right]

   &
   =
   \sum_{l = 0}^{N_y - 1}
   \vat{P}{l}
   \exp \left( \frac{2 \pi}{N_y} I \pjc l \right)
   \exp \left( + \frac{2 \pi}{N_y} I l \right).

Thus I find the left-hand side of the discrete Poisson equation consists of

.. math::

   \frac{
      \frac{
         \vat{p}{\pipp, \pjc}
         -
         \vat{p}{\pic,  \pjc}
      }{
         \vat{\Delta x}{\pip}
      }
      -
      \frac{
         \vat{p}{\pic,  \pjc}
         -
         \vat{p}{\pimm, \pjc}
      }{
         \vat{\Delta x}{\pim}
      }
   }{
      \vat{\Delta x}{\pic}
   }
   &
   =
   \frac{1}{\Delta x_{\pic} \Delta x_{\pim}}
   \sum_{l = 0}^{N_y - 1}
   \vat{P}{\pimm, l}
   \exp \left( \frac{2 \pi}{N_y} I \pjc l \right)

   &
   -
   \left(
      \frac{1}{\Delta x_{\pic} \Delta x_{\pim}}
      +
      \frac{1}{\Delta x_{\pic} \Delta x_{\pip}}
   \right)
   \sum_{l = 0}^{N_y - 1}
   \vat{P}{\pic, l}
   \exp \left( \frac{2 \pi}{N_y} I \pjc l \right)

   &
   +
   \frac{1}{\Delta x_{\pip} \Delta x_{\pip}}
   \sum_{l = 0}^{N_y - 1}
   \vat{P}{\pipp, l}
   \exp \left( \frac{2 \pi}{N_y} I \pjc l \right)

and

.. math::

   \frac{1}{\Delta y^2}
   \left(
      p_{\pic, \pjpp}
      -
      2 p_{\pic, \pjc}
      +
      p_{\pic, \pjmm}
   \right)
   &
   =
   \frac{1}{\Delta y^2}
   \sum_{l = 0}^{N_y - 1}
   \vat{P}{\pic, l}
   \exp \left( \frac{2 \pi}{N_y} I \pjc l \right)
   \left[
      \exp \left( + \frac{2 \pi}{N_y} I l \right)
      -
      2
      +
      \exp \left( - \frac{2 \pi}{N_y} I l \right)
   \right]

   &
   =
   \frac{1}{\Delta y^2}
   \sum_{l = 0}^{N_y - 1}
   \vat{P}{\pic, l}
   \exp \left( \frac{2 \pi}{N_y} I \pjc l \right)
   \left[
      2
      \cos \left( \frac{2 \pi}{N_y} l \right)
      -
      2
   \right]

   &
   =
   -
   \sum_{l = 0}^{N_y - 1}
   \vat{P}{\pic, l}
   \exp \left( \frac{2 \pi}{N_y} I \pjc l \right)
   \frac{
      4 \sin^2 \left( \frac{\pi}{N_y} l \right)
   }{
      \Delta y^2
   },

while the right-hand side yields

.. math::

   q_{\pic, \pjc}
   =
   \sum_{l = 0}^{N_y - 1}
   \vat{Q}{\pic, l}
   \exp \left( \frac{2 \pi}{N_y} I \pjc l \right).

Finally, by adopting the orthogonality of the trigonometric functions, I obtain

.. math::

   \frac{1}{\Delta x_{\pic} \Delta x_{\pim}}
   \vat{P}{\pimm, l}
   -
   \left(
      \frac{1}{\Delta x_{\pic} \Delta x_{\pim}}
      +
      \frac{1}{\Delta x_{\pic} \Delta x_{\pip}}
      +
      \frac{
         4 \sin^2 \left( \frac{\pi}{N_y} l \right)
      }{
         \Delta y^2
      }
   \right)
   \vat{P}{\pic, l}
   +
   \frac{1}{\Delta x_{\pip} \Delta x_{\pip}}
   \vat{P}{\pipp, l}
   =
   \vat{Q}{\pic, l},

where

.. math::

   P_{\pic, l}
   \equiv
   \sum_{l = 0}^{N_y - 1}
   \vat{p}{\pic, \pjc}
   \exp \left( - \frac{2 \pi}{N_y} I \pjc l \right),

   P_{\pic, l}
   \equiv
   \sum_{l = 0}^{N_y - 1}
   \vat{p}{\pic, \pjc}
   \exp \left( - \frac{2 \pi}{N_y} I \pjc l \right),

and the conversions between the physical space to the spectral space can be done using the Fast Fourier Transforms.

Note that the :math:`j` variance which was present in the original equation is now gone, and instead I obtain a linear system in :math:`x` direction, which is tri-diagonal:

.. math::

   \vat{l}{\pic} \vat{P}{\pimm, l}
   +
   \vat{c}{\pic} \vat{P}{\pic,  l}
   +
   \vat{u}{\pic} \vat{P}{\pipp, l}
   =
   \vat{Q}{\pic, l}.

See above for the definitions of the coefficients.

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

Note that the resulting total computational cost is roughly

.. math::

   \mathcal{O} \left( N_x N_y \log N_y \right),

which is a drastic reduction from the original value:

.. math::

   \mathcal{O} \left( N_x^2 N_y^2 \right).

