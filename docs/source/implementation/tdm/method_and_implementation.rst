#########################
Method and implementation
#########################

****************
Thomas algorithm
****************

========
Overview
========

Let us consider a typical tri-diagonal system whose size is :math:`n`:

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

I first look at how this matrix is solved.
For notational simplicity, I concatenate the tri-diagonal matrix and the right-hand-side term:

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

Our objective is to convert the tri-diagonal matrix to an identity matrix (i.e. matrix inversion):

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

To solve this system, there are mainly step steps: forward substitution, and the backward substitution, which are explained below.

The above system assumes Dirichlet or Neumann boundary conditions.
Extension to the periodic boundaries is also discussed below.

.. toctree::
   :maxdepth: 1

   forward
   backward
   sherman_morrison

.. note::

   Although the function implemented in this project behaves similarly as the function implemented in `LAPACK <https://netlib.org/lapack/>`_, there are mainly two differences:

   * No pivoting

      Although functions in ``LAPACK`` include pivoting operations to stabilise the system, functions implemented in this source do not since all systems to be solved in this project are (semi-) diagonally dominant.

   * Treatment of the singularity

      When a singularity is detected, functions in ``LAPACK`` terminate.
      This is not preferable in this project, since singular systems appear because of the Neumann or the periodic boundary conditions.

