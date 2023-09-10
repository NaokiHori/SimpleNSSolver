I need

   * :math:`\ux`: ``i = [1: isize + 1]`` and ``j = [0: jsize + 1]``

   * others: ``i = [0: isize + 1]`` and ``j = [0: jsize + 1]``

at ``k = 0`` and ``k = ksize + 1``.

Since they are two-dimensional (:math:`xy`) slices which are contiguous in memory, I define ``MPI_Datatype`` using ``MPI_Type_contiguous``:

.. myliteralinclude:: /../../src/halo.c
   :language: c
   :tag: define datatype in z

.. note::

   ``UZ(i, 0, ksize + 1)`` are needed to describe the advection of :math:`\uy` in the :math:`z` direction.
   Also ``UY(i, jsize + 1, 0)`` are needed to describe the advection of :math:`\uz` in the :math:`y` direction.
   Although these information are not contained by the neighbouring processes originally, they are obtained after :math:`y` halo cells are updated.

