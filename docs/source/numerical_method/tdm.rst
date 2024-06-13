
.. _tdm:

##########################
Tri-diagonal Matrix Solver
##########################

******************
Non-Periodic Cases
******************

We consider a linear system

.. math::

    \vat{l}{i}
    \vat{p}{i - 1}
    +
    \vat{c}{i}
    \vat{p}{i}
    +
    \vat{u}{i}
    \vat{p}{i + 1}
    =
    \vat{q}{i},

which appears as a consequence of the Laplace operators in the wall-normal direction.
This linear system is written as a tri-diagonal matrix:

.. math::

    \begin{bmatrix}
       c_0    & u_0    & 0      & \cdots & 0       & 0       & 0       & \cdots & 0       & 0       & 0       \\
       l_1    & c_1    & u_1    & \cdots & 0       & 0       & 0       & \cdots & 0       & 0       & 0       \\
       0      & l_2    & c_2    & \cdots & 0       & 0       & 0       & \cdots & 0       & 0       & 0       \\
       \vdots & \vdots & \vdots & \ddots & \vdots  & \vdots  & \vdots  &        & \vdots  & \vdots  & \vdots  \\
       0      & 0      & 0      & \cdots & c_{i-1} & u_{i-1} & 0       & \cdots & 0       & 0       & 0       \\
       0      & 0      & 0      & \cdots & l_{i  } & c_{i  } & u_{i  } & \cdots & 0       & 0       & 0       \\
       0      & 0      & 0      & \cdots & 0       & l_{i+1} & c_{i+1} & \cdots & 0       & 0       & 0       \\
       \vdots & \vdots & \vdots &        & \vdots  & \vdots  & \vdots  & \ddots & \vdots  & \vdots  & \vdots  \\
       0      & 0      & 0      & \cdots & 0       & 0       & 0       & \cdots & c_{n-3} & u_{n-3} & 0       \\
       0      & 0      & 0      & \cdots & 0       & 0       & 0       & \cdots & l_{n-2} & c_{n-2} & u_{n-2} \\
       0      & 0      & 0      & \cdots & 0       & 0       & 0       & \cdots & 0       & l_{n-1} & c_{n-1}
    \end{bmatrix}
    \begin{bmatrix}
       x_0     \\
       x_1     \\
       x_2     \\
       \vdots  \\
       x_{i-1} \\
       x_{i  } \\
       x_{i+1} \\
       \vdots  \\
       x_{n-3} \\
       x_{n-2} \\
       x_{n-1}
    \end{bmatrix}
    =
    \begin{bmatrix}
       q_0     \\
       q_1     \\
       q_2     \\
       \vdots  \\
       q_{i-1} \\
       q_{i  } \\
       q_{i+1} \\
       \vdots  \\
       q_{n-3} \\
       q_{n-2} \\
       q_{n-1}
    \end{bmatrix},

where :math:`x` is the answer of the system and to be computed.
In the code, a buffer ``q`` is used to store the solution :math:`x` as well as to store the input (right-hand-side terms).
Thus the input array is overwritten by the solver.

.. note::

   Although ``l[0]`` and ``u[n-1]`` are not used in this case, ``l``, ``c``, and ``u`` all have the length :math:`n` for simplicity.
   In particular these values are used for periodic systems.

We first look at how this matrix is solved.
For notational simplicity, we concatenate the tri-diagonal matrix and the right-hand-side vector:

.. math::

   \begin{bmatrix}
      c_0    & u_0    & 0      & \cdots & 0       & 0       & 0       & \cdots & 0       & 0       & 0       & q_0     \\
      l_1    & c_1    & u_1    & \cdots & 0       & 0       & 0       & \cdots & 0       & 0       & 0       & q_1     \\
      0      & l_2    & c_2    & \cdots & 0       & 0       & 0       & \cdots & 0       & 0       & 0       & q_2     \\
      \vdots & \vdots & \vdots & \ddots & \vdots  & \vdots  & \vdots  &        & \vdots  & \vdots  & \vdots  & \vdots  \\
      0      & 0      & 0      & \cdots & c_{i-1} & u_{i-1} & 0       & \cdots & 0       & 0       & 0       & q_{i-1} \\
      0      & 0      & 0      & \cdots & l_{i  } & c_{i  } & u_{i  } & \cdots & 0       & 0       & 0       & q_{i  } \\
      0      & 0      & 0      & \cdots & 0       & l_{i+1} & c_{i+1} & \cdots & 0       & 0       & 0       & q_{i+1} \\
      \vdots & \vdots & \vdots &        & \vdots  & \vdots  & \vdots  & \ddots & \vdots  & \vdots  & \vdots  & \vdots  \\
      0      & 0      & 0      & \cdots & 0       & 0       & 0       & \cdots & c_{n-3} & u_{n-3} & 0       & q_{n-3} \\
      0      & 0      & 0      & \cdots & 0       & 0       & 0       & \cdots & l_{n-2} & c_{n-2} & u_{n-2} & q_{n-2} \\
      0      & 0      & 0      & \cdots & 0       & 0       & 0       & \cdots & 0       & l_{n-1} & c_{n-1} & q_{n-1}
   \end{bmatrix},

and consider `the Gaussian elimination <https://en.wikipedia.org/wiki/Gaussian_elimination>`_.

Our objective is to convert the tri-diagonal matrix to an identity matrix (i.e., matrix inversion):

.. math::

   \begin{bmatrix}
      1      & 0      & 0      & \cdots & 0       & 0       & 0       & \cdots & 0       & 0       & 0       & x_0     \\
      0      & 1      & 0      & \cdots & 0       & 0       & 0       & \cdots & 0       & 0       & 0       & x_1     \\
      0      & 0      & 1      & \cdots & 0       & 0       & 0       & \cdots & 0       & 0       & 0       & x_2     \\
      \vdots & \vdots & \vdots & \ddots & \vdots  & \vdots  & \vdots  &        & \vdots  & \vdots  & \vdots  & \vdots  \\
      0      & 0      & 0      & \cdots & 1       & 0       & 0       & \cdots & 0       & 0       & 0       & x_{i-1} \\
      0      & 0      & 0      & \cdots & 0       & 1       & 0       & \cdots & 0       & 0       & 0       & x_{i  } \\
      0      & 0      & 0      & \cdots & 0       & 0       & 1       & \cdots & 0       & 0       & 0       & x_{i+1} \\
      \vdots & \vdots & \vdots &        & \vdots  & \vdots  & \vdots  & \ddots & \vdots  & \vdots  & \vdots  & \vdots  \\
      0      & 0      & 0      & \cdots & 0       & 0       & 0       & \cdots & 1       & 0       & 0       & x_{n-3} \\
      0      & 0      & 0      & \cdots & 0       & 0       & 0       & \cdots & 0       & 1       & 0       & x_{n-2} \\
      0      & 0      & 0      & \cdots & 0       & 0       & 0       & \cdots & 0       & 0       & 1       & x_{n-1}
   \end{bmatrix}.

=============
Forward Sweep
=============

First we try to eliminate the lower-diagonal components :math:`l_1, l_2, \cdots, l_{n-2}, l_{n-1}`, which is called the forward sweep.

---------
First row
---------

To get started, we look at the top two rows:

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

In the code, we have

.. myliteralinclude:: /../../src/tdm.c
   :language: c
   :tag: divide the first row by center-diagonal term

.. note::

   One may notice that ``u`` is modified in the above equation.
   If we implement in this way, we need to re-initialise :math:`u_i` for each right-hand-side term, which is redundant.
   To avoid this, we use a buffer ``v``, in which the modified ``u`` is stored instead of overwriting ``u``.

Next, to eliminate :math:`l_1`, we subtract *the first row times* :math:`l_1` from *the second row*, yielding

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

-------------
General Cases
-------------

We consider to extend the above process to a general :math:`i-1`-th and :math:`i`-th rows:

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
   :tag: forward sweep

--------
Last Row
--------

Basically we can do the same thing.
For the last row :math:`i = n-1`, however, the denominator

.. math::

   c_{n-1} - l_{n-1} v_{n-2}

can be :math:`0`, namely the rank of the matrix is :math:`n-1`.

This is expected, since we often impose the Neumann or the periodic boundary conditions, which can solve the differential equations only up to a constant.

In order to take into account the singularity and to avoid the resulting zero divisions, we need a special treatment:

.. myliteralinclude:: /../../src/tdm.c
   :language: c
   :tag: last row, do the same thing but consider singularity

=====================
Backward Substitution
=====================

After the forward sweep, we are left with the following system:

.. math::

   \begin{bmatrix}
      1      & v_0    & 0      & \cdots & 0       & 0       & 0      & \cdots & 0       & 0       & 0       \\
      0      & 1      & v_1    & \cdots & 0       & 0       & 0      & \cdots & 0       & 0       & 0       \\
      0      & 0      & 1      & \cdots & 0       & 0       & 0      & \cdots & 0       & 0       & 0       \\
      \vdots & \vdots & \vdots & \ddots & \vdots  & \vdots  & \vdots &        & \vdots  & \vdots  & \vdots  \\
      0      & 0      & 0      & \cdots & 1       & v_{i-1} & 0      & \cdots & 0       & 0       & 0       \\
      0      & 0      & 0      & \cdots & 0       & 1       & v_i    & \cdots & 0       & 0       & 0       \\
      0      & 0      & 0      & \cdots & 0       & 0       & 1      & \cdots & 0       & 0       & 0       \\
      \vdots & \vdots & \vdots &        & \vdots  & \vdots  & \vdots & \ddots & \vdots  & \vdots  & \vdots  \\
      0      & 0      & 0      & \cdots & 0       & 0       & 0      & \cdots & 1       & v_{n-3} & 0       \\
      0      & 0      & 0      & \cdots & 0       & 0       & 0      & \cdots & 0       & 1       & v_{n-2} \\
      0      & 0      & 0      & \cdots & 0       & 0       & 0      & \cdots & 0       & 0       & 1
   \end{bmatrix}
   \begin{bmatrix}
      x_0     \\
      x_1     \\
      x_2     \\
      \vdots  \\
      x_{i-1} \\
      x_{i  } \\
      x_{i+1} \\
      \vdots  \\
      x_{n-3} \\
      x_{n-2} \\
      x_{n-1}
   \end{bmatrix}
   =
   \begin{bmatrix}
      r_0     \\
      r_1     \\
      r_2     \\
      \vdots  \\
      r_{i-1} \\
      r_{i  } \\
      r_{i+1} \\
      \vdots  \\
      r_{n-3} \\
      r_{n-2} \\
      r_{n-1}
   \end{bmatrix}.

In the last row, we have

.. math::

   x_{n-1} = r_{n-1},

which has already been computed in the forward sweep.

Also, since we have

.. math::

   x_i = r_i - v_i x_{i+1},

we can compute :math:`x_i` one after another (sequentially from :math:`i = n-2` to :math:`i = 0`):

.. myliteralinclude:: /../../src/tdm.c
   :language: c
   :tag: backward substitution

Note again that ``q`` is shared among :math:`x_i` (output) and :math:`q_i` (input) in the code.

**************
Periodic Cases
**************

The original Thomas algorithm only considers the tri-diagonal matrix.
With periodic boundary conditions, right-top and left-bottom corners have non-zero values (see below).

Fortunately, by using `the Sherman-Morrison formula <https://en.wikipedia.org/wiki/Sherman–Morrison_formula>`_, we can handle this minor correction in the framework of the Thomas algorithm.

Now, we consider the following system:

.. math::

    \begin{bmatrix}
       c_0     & u_0    & 0      & \cdots & 0       & 0       & 0       & \cdots & 0       & 0       & l_0     \\
       l_1     & c_1    & u_1    & \cdots & 0       & 0       & 0       & \cdots & 0       & 0       & 0       \\
       0       & l_2    & c_2    & \cdots & 0       & 0       & 0       & \cdots & 0       & 0       & 0       \\
       \vdots  & \vdots & \vdots & \ddots & \vdots  & \vdots  & \vdots  &        & \vdots  & \vdots  & \vdots  \\
       0       & 0      & 0      & \cdots & c_{i-1} & u_{i-1} & 0       & \cdots & 0       & 0       & 0       \\
       0       & 0      & 0      & \cdots & l_{i  } & c_{i  } & u_{i  } & \cdots & 0       & 0       & 0       \\
       0       & 0      & 0      & \cdots & 0       & l_{i+1} & c_{i+1} & \cdots & 0       & 0       & 0       \\
       \vdots  & \vdots & \vdots &        & \vdots  & \vdots  & \vdots  & \ddots & \vdots  & \vdots  & \vdots  \\
       0       & 0      & 0      & \cdots & 0       & 0       & 0       & \cdots & c_{n-3} & u_{n-3} & 0       \\
       0       & 0      & 0      & \cdots & 0       & 0       & 0       & \cdots & l_{n-2} & c_{n-2} & u_{n-2} \\
       u_{n-1} & 0      & 0      & \cdots & 0       & 0       & 0       & \cdots & 0       & l_{n-1} & c_{n-1}
    \end{bmatrix}
    \begin{bmatrix}
       x_0     \\
       x_1     \\
       x_2     \\
       \vdots  \\
       x_{i-1} \\
       x_{i  } \\
       x_{i+1} \\
       \vdots  \\
       x_{n-3} \\
       x_{n-2} \\
       x_{n-1}
    \end{bmatrix}
    =
    \begin{bmatrix}
       q_0     \\
       q_1     \\
       q_2     \\
       \vdots  \\
       q_{i-1} \\
       q_{i  } \\
       q_{i+1} \\
       \vdots  \\
       q_{n-3} \\
       q_{n-2} \\
       q_{n-1}
    \end{bmatrix},

where one may notice that the top-right and the bottom-left corners have non-zero values.

Since we have

.. math::

   \begin{alignat}{5}
      & c_0     x_0     & & + u_0     x_1     & & + l_0     x_{n-1} & & = q_0     & \,\,\, &        0           \text{-th row} \\
      & l_{n-2} x_{n-3} & & + c_{n-2} x_{n-2} & & + u_{n-2} x_{n-1} & & = q_{n-2} & \,\,\, & \left( n-2 \right) \text{-th row}
   \end{alignat}

or

.. math::

   \begin{alignat}{3}
      & c_0     x_0     & & + u_0     x_1     & & = q_0     - l_0     x_{n-1} \\
      & l_{n-2} x_{n-3} & & + c_{n-2} x_{n-2} & & = q_{n-2} - u_{n-2} x_{n-1},
   \end{alignat}

we can *shrink* the system:

.. math::

   \begin{bmatrix}
      c_0     & u_0    & 0      & \cdots & 0       & 0       & 0       & \cdots & 0       & 0       & 0       \\
      l_1     & c_1    & u_1    & \cdots & 0       & 0       & 0       & \cdots & 0       & 0       & 0       \\
      0       & l_2    & c_2    & \cdots & 0       & 0       & 0       & \cdots & 0       & 0       & 0       \\
      \vdots  & \vdots & \vdots & \ddots & \vdots  & \vdots  & \vdots  &        & \vdots  & \vdots  & \vdots  \\
      0       & 0      & 0      & \cdots & c_{i-1} & u_{i-1} & 0       & \cdots & 0       & 0       & 0       \\
      0       & 0      & 0      & \cdots & l_{i  } & c_{i  } & u_{i  } & \cdots & 0       & 0       & 0       \\
      0       & 0      & 0      & \cdots & 0       & l_{i+1} & c_{i+1} & \cdots & 0       & 0       & 0       \\
      \vdots  & \vdots & \vdots &        & \vdots  & \vdots  & \vdots  & \ddots & \vdots  & \vdots  & \vdots  \\
      0       & 0      & 0      & \cdots & 0       & 0       & 0       & \cdots & c_{n-4} & u_{n-4} & 0       \\
      0       & 0      & 0      & \cdots & 0       & 0       & 0       & \cdots & l_{n-3} & c_{n-3} & u_{n-3} \\
      0       & 0      & 0      & \cdots & 0       & 0       & 0       & \cdots & 0       & l_{n-2} & c_{n-2}
   \end{bmatrix}
   \begin{bmatrix}
      x_0     \\
      x_1     \\
      x_2     \\
      \vdots  \\
      x_{i-1} \\
      x_{i  } \\
      x_{i+1} \\
      \vdots  \\
      x_{n-4} \\
      x_{n-3} \\
      x_{n-2}
   \end{bmatrix}
   =
   \begin{bmatrix}
      q_0 - l_0 x_{n-1} \\
      q_1               \\
      q_2               \\
      \vdots            \\
      q_{i-1}           \\
      q_{i  }           \\
      q_{i+1}           \\
      \vdots            \\
      q_{n-4}           \\
      q_{n-3}           \\
      q_{n-2} - u_{n-2} x_{n-1}
   \end{bmatrix},

i.e., we moved :math:`x_{n-1}` from the left-hand side to the right-hand side, and the size of the system is now :math:`n - 1`.

For notational simplicity, hereafter we write this as

.. math::

   \underline{\underline{A}} \, \underline{x} = \underline{q}.

One may notice that this treatment has removed the additional components in the original matrix coming from the periodicity, and as a result we go back to the tri-diagonal system.

The new system, however, includes an unknown :math:`x_{n-1}`.
As soon as we try to start the forward sweep, we would be in trouble since the first row in the right-hand side includes unknown value.
To resolve this situation, we consider to split the system into two problems:

.. math::

   {\underline{\underline{A}}} \, {\underline{x}}^0 & = {\underline{q}}^0, \\
   {\underline{\underline{A}}} \, {\underline{x}}^1 & = {\underline{q}}^1,

where we define

.. math::

   \underline{q}^0
   =
   \begin{bmatrix}
      q_0     \\
      q_1     \\
      q_2     \\
      \vdots  \\
      q_{i-1} \\
      q_{i  } \\
      q_{i+1} \\
      \vdots  \\
      q_{n-4} \\
      q_{n-3} \\
      q_{n-2}
   \end{bmatrix},
   \underline{q}^1
   =
   \begin{bmatrix}
      - l_0  \\
      0      \\
      0      \\
      \vdots \\
      0      \\
      0      \\
      0      \\
      \vdots \\
      0      \\
      0      \\
      - u_{n-2}
   \end{bmatrix},

which satisfies

.. math::

   {\underline{q}}
   =
   {\underline{q}}^0
   +
   x_{n-1}
   {\underline{q}}^1.

Note that the superscripts are used to distinguish the two problems (not the exponents).

Since these two systems:

.. math::

   {\underline{\underline{A}}} \, {\underline{x}}^0 & = {\underline{q}}^0, \\
   {\underline{\underline{A}}} \, {\underline{x}}^1 & = {\underline{q}}^1,

do not contain any unknown, we can solve them as two independent tri-diagonal systems:

.. math::

   {\underline{x}}
   =
   {\underline{x}}^0
   +
   x_{n-1}
   \times
   {\underline{x}}^1,

indicating that, the solution of the original system is the superposition of the solutions of the two tri-diagonal systems.

The last piece is how to find :math:`x_{n-1}`, which is obtained by looking at the relation:

.. math::

   u_{n-1} x_0 + l_{n-1} x_{n-2} + c_{n-1} x_{n-1} = q_{n-1},

which appears in the last row of the original (:math:`n \times n`) system.

Since we have

.. math::

   x_0     & = x_0^0     + x_{n-1} \times x_0^1,     \\
   x_{n-2} & = x_{n-2}^0 + x_{n-1} \times x_{n-2}^1, \\

we notice

.. math::

     \left( u_{n-1} x_0^0 + l_{n-1} x_{n-2}^0           \right)
   + \left( u_{n-1} x_0^1 + l_{n-1} x_{n-2}^1 + c_{n-1} \right) x_{n-1}
   = q_{n-1},

and thus

.. math::

   x_{n-1} =
      \frac{q_{n-1} - u_{n-1} x_0^0 - l_{n-1} x_{n-2}^0}
      {c_{n-1} + u_{n-1} x_0^1 + l_{n-1} x_{n-2}^1}.

This relation indicates that :math:`x_{n-1}` can be computed after solving the two shrunk linear systems

.. math::

   {\underline{\underline{A}}} \, {\underline{x}}^0 = {\underline{q}}^0`,

.. math::

   {\underline{\underline{A}}} \, {\underline{x}}^1 = {\underline{q}}^1`.

Here is the summary and the corresponding implementation:

#. Solve :math:`{\underline{\underline{A}}} \, {\underline{x}}^1 = {\underline{q}}^1`:

   .. myliteralinclude:: /../../src/tdm.c
      :language: c
      :tag: solve additional system coming from periodicity

#. Solve :math:`{\underline{\underline{A}}} \, {\underline{x}}^0 = {\underline{q}}^0`:

   .. myliteralinclude:: /../../src/tdm.c
      :language: text
      :tag: solve normal system

   .. note::

      The input argument ``q`` includes multiple (``nrhs``, corresponding to the loop whose index is ``j``) right-hand-side terms.
      We assign the pointer of each right-hand side to ``q0`` here.

#. Find :math:`x_{n-1}`

   :math:`x_{n-1}` is updated following

   .. math::

      x_{n-1} = \frac{q_{n-1} - u_{n-1} x_0^0 - l_{n-1} x_{n-2}^0}{c_{n-1} + u_{n-1} x_0^1 + l_{n-1} x_{n-2}^1}:

   .. myliteralinclude:: /../../src/tdm.c
      :language: c
      :tag: find x_{n-1}

#. Compute the solution of the original system :math:`{\underline{\underline{A}}} \, {\underline{x}} = {\underline{q}}`

   We use

   .. math::

      {\underline{x}} = {\underline{x}}^0 + x_{n-1} \times {\underline{x}}^1:

   .. myliteralinclude:: /../../src/tdm.c
      :language: c
      :tag: solve original system

.. note::

   Although the function implemented in this project behaves similarly as the function implemented in `LAPACK <https://netlib.org/lapack/>`_, there are mainly two differences:

   * No pivoting

      Although functions in ``LAPACK`` include pivoting operations to stabilise the system, functions implemented in this source do not since all systems to be solved in this project are (semi-) diagonally dominant.

   * Treatment of the singularity

      When a singularity is detected, functions in ``LAPACK`` terminate.
      This is not preferable in this project, since singular systems appear because of the Neumann or the periodic boundary conditions.

