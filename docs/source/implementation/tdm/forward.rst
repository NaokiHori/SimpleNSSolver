
.. _forward_substitution:

####################
Forward substitution
####################

First I try to eliminate the lower-diagonal components :math:`l_1, l_2, \cdots, l_{n-2}, l_{n-1}`, which is called the forward substitution.

****************
:math:`i = 0, 1`
****************

To get started, I look at the top two rows:

.. math::

   \begin{bmatrix}
      c_0 & u_0 & 0   & \cdots & 0 & 0 & 0 & \cdots & 0 & 0 & 0 & q_0 \\
      l_1 & c_1 & u_1 & \cdots & 0 & 0 & 0 & \cdots & 0 & 0 & 0 & q_1
   \end{bmatrix}.

Dividing the first row by :math:`c_0` yields

.. math::

   &
   \begin{bmatrix}
      1   & u_0 / c_0 & 0   & \cdots & 0 & 0 & 0 & \cdots & 0 & 0 & 0 & q_0 / c_0 \\
      l_1 & c_1       & u_1 & \cdots & 0 & 0 & 0 & \cdots & 0 & 0 & 0 & q_1
   \end{bmatrix} \\
   &
   = \\
   &
   \begin{bmatrix}
      1   & v_0 & 0   & \cdots & 0 & 0 & 0 & \cdots & 0 & 0 & 0 & r_0 \\
      l_1 & c_1 & u_1 & \cdots & 0 & 0 & 0 & \cdots & 0 & 0 & 0 & q_1
   \end{bmatrix},

namely

.. math::

   v_0 \leftarrow \frac{u_0}{c_0},

.. math::

   r_0 \leftarrow \frac{q_0}{c_0}.

In the code, I have

.. myliteralinclude:: /../../src/tdm.c
   :language: c
   :tag: divide the first row by center-diagonal term

.. note::

   One may notice that ``u`` is modified in the above equation.
   If I implement in this way, I need to re-initialise :math:`u_i` for each right-hand-side term, which is redundant.
   To avoid this, I use a buffer ``v``, in which the modified ``u`` is stored instead of overwriting ``u``.

Next, to eliminate :math:`l_1`, I subtract *the first row times* :math:`l_1` from *the second row*, yielding

.. math::

   \begin{bmatrix}
      1 & v_0           & 0   & \cdots & 0 & 0 & 0 & \cdots & 0 & 0 & 0 & r_0 \\
      0 & c_1 - l_1 v_0 & u_1 & \cdots & 0 & 0 & 0 & \cdots & 0 & 0 & 0 & q_1 - l_1 r_0
   \end{bmatrix},

or

.. math::

   &
   \begin{bmatrix}
      1 & v_0 & 0                                  & \cdots & 0 & 0 & 0 & \cdots & 0 & 0 & 0 & r_0 \\
      0 & 1            & \frac{u_1}{c_1 - l_1 v_0} & \cdots & 0 & 0 & 0 & \cdots & 0 & 0 & 0 & \frac{q_1 - l_1 r_0}{c_1 - l_1 v_0}
   \end{bmatrix}, \\
   &
   = \\
   &
   \begin{bmatrix}
      1 & v_0 & 0   & \cdots & 0 & 0 & 0 & \cdots & 0 & 0 & 0 & r_0 \\
      0 & 1   & v_1 & \cdots & 0 & 0 & 0 & \cdots & 0 & 0 & 0 & r_1
   \end{bmatrix}.

***********************************
General :math:`(i-1)` and :math:`i`
***********************************

I consider to extend the above process to a general :math:`i-1`-th and :math:`i`-th rows:

.. math::

   \begin{bmatrix}
      0 & 0 & 0 & \cdots & 1   & v_{i-1} & 0   & \cdots & 0 & 0 & 0 & r_{i-1} \\
      0 & 0 & 0 & \cdots & l_i & c_i     & u_i & \cdots & 0 & 0 & 0 & q_{i  }
   \end{bmatrix},

where the upper row (:math:`i-1`-th row) has already been updated, while the bottom row (:math:`i`-th row) is to be updated now.

Now let us consider to eliminate :math:`l_i`:

.. math::

   \begin{bmatrix}
      0 & 0 & 0 & \cdots & 1   & v_{i-1}           & 0   & \cdots & 0 & 0 & 0 & r_{i-1} \\
      0 & 0 & 0 & \cdots & 0   & c_i - l_i v_{i-1} & u_i & \cdots & 0 & 0 & 0 & q_{i} - l_i r_{i-1}
   \end{bmatrix},

or

.. math::

   &
   \begin{bmatrix}
      0 & 0 & 0 & \cdots & 1 & v_{i-1} & 0                             & \cdots & 0 & 0 & 0 & r_{i-1} \\
      0 & 0 & 0 & \cdots & 0 & 1       & \frac{u_i}{c_i - l_i v_{i-1}} & \cdots & 0 & 0 & 0 & \frac{q_{i} - l_i r_{i-1}}{c_i - l_i v_{i-1}}
   \end{bmatrix} \\
   &
   = \\
   &
   \begin{bmatrix}
      0 & 0 & 0 & \cdots & 1   & v_{i-1} & 0     & \cdots & 0 & 0 & 0 & r_{i-1} \\
      0 & 0 & 0 & \cdots & 0   & 1       & v_{i} & \cdots & 0 & 0 & 0 & r_{i  }
   \end{bmatrix},

namely,

.. math::

   v_i \leftarrow \frac{u_i}{c_i - l_i v_{i-1}},

.. math::

   r_i \leftarrow \frac{q_i - l_i r_{i-1}}{c_i - l_i v_{i-1}}.

This is repeated from :math:`i = 1` to :math:`n - 2`:

.. myliteralinclude:: /../../src/tdm.c
   :language: c
   :tag: forward substitution

*****************
:math:`i = n - 1`
*****************

Basically I can do the same thing.
For the last row :math:`i = n-1`, however, the denominator

.. math::

   c_{n-1} - l_{n-1} v_{n-2}

can be :math:`0`, namely the rank of the matrix is :math:`n-1`.

This is expected, since I often impose the Neumann or the periodic boundary conditions, which can only solve the differential equations up to a constant.

In order to take into account the singularity and to avoid the resulting zero divisions, I need a special treatment:

.. myliteralinclude:: /../../src/tdm.c
   :language: c
   :tag: last row, do the same thing but consider singularity

.. seealso::

   This forward substitution is followed by :ref:`the backward substitution <backward_substitution>`.

