
.. _poisson_eigen_decomposition:

###################
Eigen-decomposition
###################

*************************
Discrete Poisson equation
*************************

By adopting the second-order-accurate central-difference scheme, the two-dimensional Poisson equation formulated in :ref:`the previous section <poisson_governing_equation>` is discretised as

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

Note that I assume the grid sizes in the :math:`y` direction is uniform, while those in the :math:`x` direction can vary.

This yields a linear system (matrix) whose size is :math:`\left( N_x \times N_y \right)^2`, where :math:`N_x` and :math:`N_y` are the number of cell centers in each direction.
Although the resulting matrix is sparse, it is still challenging to solve it by means of a direct or an iterative method.

*******************
Eigen-decomposition
*******************

To solve the system efficiently, I consider the eigen decomposition of the system in the uniform direction (:math:`y`), i.e. I aim to reduce the scattered nature of the original linear system by diagonalising it in the :math:`y` direction.

Since the grid points are equidistant, I choose the trigonometric functions as bases and consider the discrete Fourier series expansion defined as

.. math::

   \mathcal{F}_d \left[ \vat{p}{\pjc} \right]
   \equiv
   \sum_{j = 0}^{N_y - 1} \vat{p}{\pjc} \exp \left( - \frac{2 \pi}{N_y} I k_y \pjc \right),

where :math:`I` is used to denote the imaginary unit, while :math:`k_y` is the (discrete) wave numbers

.. math::

   k_y = 0, 1, 2, \cdots, \frac{N_y}{2} - 1,

since I assume that :math:`p \in \mathbb{R}`.

Note that, for notational convenience, I drop the indices in the :math:`x` direction in this part.

To know the relation between :math:`\vat{p}{\pjc}` and :math:`\vat{p}{\pjpp}`, I consider

.. math::

   \sum_{j = 0}^{N_y - 1} \vat{p}{\pjpp} \exp \left( - \frac{2 \pi}{N_y} I k_y \pjc \right).

By replacing :math:`j` with :math:`j - 1`, I find that the above equation is equivalent to

.. math::

   &
   \sum_{j - 1 = 0}^{N_y - 1} \vat{p}{\pjc} \exp \left\{ - \frac{2 \pi}{N_y} I k_y \left( \pjmm \right) \right\} \\
   =
   &
   \sum_{j = 1}^{N_y} \vat{p}{\pjc} \exp \left( - \frac{2 \pi}{N_y} I k_y \pjc \right) \exp \left( \frac{2 \pi}{N_y} I k_y \right).

Regarding the last term of the summation, I have

.. math::

   \vat{p}{N_y}
   =
   \vat{p}{0}

and

.. math::

   \exp \left( - \frac{2 \pi}{N_y} I k_y \times N_y \right)
   =
   \exp \left( - \frac{2 \pi}{N_y} I k_y \times   0 \right)

and thus

.. math::

   \exp \left( \frac{2 \pi}{N_y} I k_y \right) \sum_{j = 0}^{N_y - 1} \vat{p}{\pjc} \exp \left( - \frac{2 \pi}{N_y} I k_y \pjc \right)
   =
   \exp \left( \frac{2 \pi}{N_y} I k_y \right) \mathcal{F}_d \left[ \vat{p}{\pjc} \right].

In summary, I have

.. math::

   \mathcal{F}_d \left[ \vat{p}{\pjpp} \right]
   =
   \exp \left( + \frac{2 \pi}{N_y} I k_y \right) \mathcal{F}_d \left[ \vat{p}{\pjc} \right],

and similarly

.. math::

   \mathcal{F}_d \left[ \vat{p}{\pjmm} \right]
   =
   \exp \left( - \frac{2 \pi}{N_y} I k_y \right) \mathcal{F}_d \left[ \vat{p}{\pjc} \right].

Finally, because of the linearity, I obtain

.. math::

   \mathcal{F}_d \left[
      \vat{p}{\pjpp}
      -
      2 \vat{p}{\pjc}
      +
      \vat{p}{\pjmm}
   \right]
   =
   \left[
      \exp \left( + \frac{2 \pi}{N_y} I k_y \right)
      +
      \exp \left( - \frac{2 \pi}{N_y} I k_y \right)
      -
      2
   \right]
   \times
   \mathcal{F}_d \left[ \vat{p}{\pjc} \right],

or thanks to the Euler's formula:

.. math::

   \left[
      2 \cos \left( \frac{2 \pi}{N_y} k_y \right)
      -
      2
   \right]
   \times
   \mathcal{F}_d \left[ \vat{p}{\pjc} \right],

or

.. math::

   \left[
      -
      4
      \sin^2 \left( \frac{\pi}{N_y} k_y \right)
   \right]
   \times
   \mathcal{F}_d \left[ \vat{p}{\pjc} \right].

As I see here, the :math:`j` variance, which exists in the original Poisson equation, is now gone, indicating that the eigen decomposition in the :math:`y` direction is completed.

Thanks to this decomposition, now I am left with a system with respect to the :math:`x` direction for each wave number :math:`k_y`, whose total size is roughly

.. math::

   N_x \times \frac{N_y}{2},

which clearly shows a drastic reduction in the computational effort, which was originally

.. math::

   \left( N_x \times N_y \right)^2

if treated naively.

In :ref:`the next section <poisson_implementation>`, the resulting system is discussed more, and the implementation is briefly described.

