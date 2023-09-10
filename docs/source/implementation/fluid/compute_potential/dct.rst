####################################
`src/fluid/compute_potential/dct.c`_
####################################

.. _src/fluid/compute_potential/dct.c: https://github.com/NaokiHori/SimpleNSSolver/blob/main/src/fluid/compute_potential/dct.c

Since I assume the spacing is uniform, I can adopt the ``FFT``-based eigen decomposition in the :math:`x` direction as well as in the other homogeneous directions.

It is possible to adopt the discrete fast Fourier transforms in all directions whose cost is :math:`\mathcal{O} \left( N \log N \right)` in each direction.
Only for the last direction, however, I can solve the equation using the Thomas algorithm whose cost is :math:`\mathcal{O} \left( N \right)`.
Although the difference should be small, I use the Thomas algorithm since the linear solver is implemented already and used in the other parts of the project.

#. Project :math:`x` direction to wave space

   .. math::

      f \left( x, y \right)
      \rightarrow
      f \left( k_x, y \right),

   .. math::

      f \left( x, y, z \right)
      \rightarrow
      f \left( k_x, y, z \right).

   .. myliteralinclude:: /../../src/fluid/compute_potential/dct.c
      :language: c
      :tag: project x to wave space

#. Transpose ``x1pencils`` to ``y1pencils``

   Memory order is changed from :math:`x`-contiguous to :math:`y`-contiguous to benefit the following procedures.

   .. myliteralinclude:: /../../src/fluid/compute_potential/dct.c
      :language: c
      :tag: transpose real x1pencil to y1pencil

#. (3D only) Project :math:`y` direction to wave space

   .. math::

      f \left( k_x, y, z \right)
      \rightarrow
      f \left( k_x, k_y, z \right).

   .. myliteralinclude:: /../../src/fluid/compute_potential/dct.c
      :language: c
      :tag: project y to wave space

#. (3D only) Transpose ``y1pencils`` to ``z1pencils``

   Memory order is changed from :math:`y`-contiguous to :math:`z`-contiguous to benefit the following procedures.

   .. myliteralinclude:: /../../src/fluid/compute_potential/dct.c
      :language: c
      :tag: transpose complex y1pencil to z1pencil

#. Solve linear systems

   I am left with the linear systems (tri-diagonal matrices), which are solved here.

   .. myliteralinclude:: /../../src/fluid/compute_potential/dct.c
      :language: c
      :tag: solve linear systems

#. (3D only) Transpose ``z1pencils`` to ``y1pencils``

   Memory order is changed from :math:`z`-contiguous to :math:`y`-contiguous to benefit the following procedures.

   .. myliteralinclude:: /../../src/fluid/compute_potential/dct.c
      :language: c
      :tag: transpose complex z1pencil to y1pencil

#. (3D only) Project :math:`y` direction to physical space

   .. math::

      f \left( k_x, k_y, z \right)
      \rightarrow
      f \left( k_x, y, z \right)

   .. myliteralinclude:: /../../src/fluid/compute_potential/dct.c
      :language: c
      :tag: project y to physical space

#. Transpose ``y1pencils`` to ``x1pencils``

   Memory order is changed from :math:`y`-contiguous to :math:`x`-contiguous to benefit the following procedures.

   .. myliteralinclude:: /../../src/fluid/compute_potential/dct.c
      :language: c
      :tag: transpose real y1pencil to x1pencil

#. Project :math:`x` direction to physical space

   .. math::

      f \left( k_x, y \right)
      \rightarrow
      f \left( x, y \right),

   .. math::

      f \left( k_x, y, z \right)
      \rightarrow
      f \left( x, y, z \right).

   .. myliteralinclude:: /../../src/fluid/compute_potential/dct.c
      :language: c
      :tag: project x to physical space

