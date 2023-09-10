
.. _sherman_morrison:

########################
Sherman-Morrison formula
########################

The original Thomas algorithm only considers the tri-diagonal matrix.
With periodic boundary conditions, right-top and left-bottom corners have non-zero values (see below).

Fortunately, by using `the Sherman-Morrison formula <https://en.wikipedia.org/wiki/Sherman–Morrison_formula>`_, I can handle this minor correction in the framework of the Thomas algorithm.

Now, I consider the following system:

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

Since I have

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

I can *shrink* the system:

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

i.e. I moved :math:`x_{n-1}` from the left-hand side to the right-hand side, and the size of the system is now :math:`n - 1`.

For notational simplicity, hereafter I write this as

.. math::

   \underline{\underline{A}} \, \underline{x} = \underline{q}.

One may notice that this treatment has removed the additional components in the original matrix coming from the periodicity, and as a result I go back to the tri-diagonal system.

The new system, however, includes an unknown :math:`x_{n-1}`.
As soon as I try to start the forward substitution, I would be in trouble since the first row in the right-hand side includes unknown value.
To resolve this situation, I consider to split the system into two problems:

.. math::

   {\underline{\underline{A}}} \, {\underline{x}}^0 & = {\underline{q}}^0, \\
   {\underline{\underline{A}}} \, {\underline{x}}^1 & = {\underline{q}}^1,

where I define

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

do not contain any unknown, I can solve them as two independent tri-diagonal systems:

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

Since I have

.. math::

   x_0     & = x_0^0     + x_{n-1} \times x_0^1,     \\
   x_{n-2} & = x_{n-2}^0 + x_{n-1} \times x_{n-2}^1, \\

I notice

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
      I assign the pointer of each right-hand side to ``q0`` here.

#. Find :math:`x_{n-1}`

   :math:`x_{n-1}` is updated following

   .. math::

      x_{n-1} = \frac{q_{n-1} - u_{n-1} x_0^0 - l_{n-1} x_{n-2}^0}{c_{n-1} + u_{n-1} x_0^1 + l_{n-1} x_{n-2}^1}:

   .. myliteralinclude:: /../../src/tdm.c
      :language: c
      :tag: find x_{n-1}

#. Compute the solution of the original system :math:`{\underline{\underline{A}}} \, {\underline{x}} = {\underline{q}}`

   I use

   .. math::

      {\underline{x}} = {\underline{x}}^0 + x_{n-1} \times {\underline{x}}^1:

   .. myliteralinclude:: /../../src/tdm.c
      :language: c
      :tag: solve original system

