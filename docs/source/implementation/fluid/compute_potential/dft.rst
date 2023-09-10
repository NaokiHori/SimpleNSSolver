####################################
`src/fluid/compute_potential/dft.c`_
####################################

.. _src/fluid/compute_potential/dft.c: https://github.com/NaokiHori/SimpleNSSolver/blob/main/src/fluid/compute_potential/dft.c

Since I assume the spacing can be non-uniform, I cannot adopt the ``FFT``-based eigen decomposition in the :math:`x` direction.
Thus, I first treat the other homogeneous directions by means of ``DFT`` to diagonalise the systems, which is followed by solving the linear systems in the :math:`x` direction.

#. Transpose ``x1pencils`` to ``y1pencils``

   .. myliteralinclude:: /../../src/fluid/compute_potential/dft.c
      :language: c
      :tag: transpose real x1pencil to y1pencil

#. Project :math:`y` direction to wave space

   .. math::

      f \left( x, y \right)
      \rightarrow
      f \left( x, k_y \right),

   .. math::

      f \left( x, y, z \right)
      \rightarrow
      f \left( x, k_y, z \right).

   .. myliteralinclude:: /../../src/fluid/compute_potential/dft.c
      :language: c
      :tag: project y to wave space

#. (2D only) Transpose ``y1pencils`` to ``x1pencils``

   .. myliteralinclude:: /../../src/fluid/compute_potential/dft.c
      :language: c
      :tag: transpose complex y1pencil to x1pencil

#. (3D only) Transpose ``y1pencils`` to ``z1pencils``

   .. myliteralinclude:: /../../src/fluid/compute_potential/dft.c
      :language: c
      :tag: transpose complex y1pencil to z1pencil

#. (3D only) Project :math:`z` direction to wave space

   .. math::

      f \left( x, k_y, z \right)
      \rightarrow
      f \left( x, k_y, k_z \right).

   .. myliteralinclude:: /../../src/fluid/compute_potential/dft.c
      :language: c
      :tag: project z to wave space

#. (3D only) Transpose ``z1pencils`` to ``x2pencils``

   .. myliteralinclude:: /../../src/fluid/compute_potential/dft.c
      :language: c
      :tag: transpose complex z1pencil to x2pencil

#. Solve linear systems

   I am left with the linear systems (tri-diagonal matrices), which are solved here.

   .. myliteralinclude:: /../../src/fluid/compute_potential/dft.c
      :language: c
      :tag: transpose complex y1pencil to z1pencil

#. (2D only) Transpose ``x1pencils`` to ``y1pencils``

   .. myliteralinclude:: /../../src/fluid/compute_potential/dft.c
      :language: c
      :tag: transpose complex x1pencil to y1pencil

#. (3D only) Transpose ``x2pencils`` to ``z1pencils``

   .. myliteralinclude:: /../../src/fluid/compute_potential/dft.c
      :language: c
      :tag: transpose complex x2pencil to z1pencil

#. (3D only) Project :math:`z` direction to physical space

   .. math::

      f \left( x, k_y, k_z \right)
      \rightarrow
      f \left( x, k_y, z \right).

   .. myliteralinclude:: /../../src/fluid/compute_potential/dft.c
      :language: c
      :tag: project z to physical space

#. (3D only) Transpose ``z1pencils`` to ``y1pencils``

   .. myliteralinclude:: /../../src/fluid/compute_potential/dft.c
      :language: c
      :tag: transpose complex z1pencil to y1pencil

#. Project :math:`y` direction to physical space

   .. math::

      f \left( x, k_y \right)
      \rightarrow
      f \left( x, y \right),

   .. math::

      f \left( x, k_y, z \right)
      \rightarrow
      f \left( x, y, z \right).

   .. myliteralinclude:: /../../src/fluid/compute_potential/dft.c
      :language: c
      :tag: project y to physical space

#. Transpose ``y1pencils`` to ``x1pencils``

   .. myliteralinclude:: /../../src/fluid/compute_potential/dft.c
      :language: c
      :tag: transpose real y1pencil to x1pencil

