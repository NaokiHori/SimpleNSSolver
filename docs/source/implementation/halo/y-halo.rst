* Two-dimensional case

   I need

      * :math:`\ux`: ``i = [1: isize + 1]``

      * others: ``i = [0: isize + 1]``

   at ``j = 0`` and ``j = jsize+1``, which are contained by the process in the negative and positive :math:`y` directions, respectively.
   Since these elements are contiguous in the memory, I define an ``MPI_Datatype`` using ``MPI_Type_contiguous``.

* Three-dimensional case

   I need

      * :math:`\ux`: ``i = [1: isize + 1]`` and ``k = [0: ksize + 1]``

      * :math:`\uy`, :math:`\uz`, :math:`p`, :math:`T`: ``i = [0: isize + 1]`` and ``k = [0: ksize + 1]``

   at ``j = 0`` and ``j = jsize + 1``.

   Since

      * there are ``ksize + 2`` blocks,

      * each block has contiguous elements (``isize + 1`` for :math:`\ux`, ``isize + 2`` otherwise),

      * stride between each block is ``(isize + 1) * (jsize + 2)`` for :math:`\ux`, ``(isize + 2) * (jsize + 2)`` otherwise,

   I define ``MPI_Datatype`` using ``MPI_Type_vector``.

.. myliteralinclude:: /../../src/halo.c
   :language: c
   :tag: define datatype in y

