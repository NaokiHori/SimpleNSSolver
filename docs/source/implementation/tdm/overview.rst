########
Overview
########

This source file contains several functions to solve tri-diagonal linear systems (I name the title based on the abbreviation of the ``tri-diagonal matrix``).
As usual, one needs to follow the three steps:

#. Initialisation

   Internal buffers are allocated.

#. Execution

   The systems are solved for specific right-hand-side terms.

#. Finalisation

   Buffers are deallocated.

**************
Initialisation
**************

.. mydeclare:: /../../src/tdm.c
   :language: c
   :tag: construct

This step is to initialise the following structure which stores all information about the linear systems to be solved:

.. myliteralinclude:: /../../src/tdm.c
   :language: c
   :tag: definition of tdm_info_t_

Once this constructor is called, one can assign the information about the tri-diagonal matrix:

.. mydeclare:: /../../src/tdm.c
   :language: c
   :tag: get_l

.. mydeclare:: /../../src/tdm.c
   :language: c
   :tag: get_c

.. mydeclare:: /../../src/tdm.c
   :language: c
   :tag: get_u

Example:

.. code-block:: c

   // call constructor
   tdm_info_t *info = null;
   int retval = tdm.construct(size, nrhs, is_periodic, is_complex, &info);
   if(0 != retval){
      printf("plan creation failed\n");
      exit(1);
   }
   // initialise tri-diagonal matrix
   double *l = null;
   double *c = null;
   double *u = null;
   tdm.get_l(info, l);
   tdm.get_c(info, c);
   tdm.get_u(info, u);
   for(int i = 0; i < size; i++){
      l[i] = set_lower_diagonals(i);
      c[i] = set_centr_diagonals(i);
      u[i] = set_upper_diagonals(i);
   }

.. note::

   ``l``, ``c`` and ``u`` are not changed in the following procedures and thus one can reuse them.

************
Solve system
************

.. mydeclare:: /../../src/tdm.c
   :language: c
   :tag: solve

By calling this function (``tdm.solve(...)``), the linear systems are solved.

Example:

.. code-block:: c

   // initialise right-hand-side term(s)
   for(int j = 0; j < nrhs; j++){
      for(int i = 0; i < size; i++){
         data[j * size + i] = set_right_hand_side(i, j);
      }
   }
   // call API
   tdm.solve(info, data);

Note that, as long as the tri-diagonal matrix is identical (``n``, ``l``, ``c`` and ``u`` are the same), one can solve multiple (``nrhs``) inputs simultaneously.

************
Finalisation
************

.. mydeclare:: /../../src/tdm.c
   :language: c
   :tag: destruct

If ``info`` is no longer used, call this destructor to free all internal memory:

.. code-block:: c

   tdm.destruct(info);

