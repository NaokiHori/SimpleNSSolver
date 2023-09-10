
.. _backward_substitution:

#####################
Backward substitution
#####################

After :ref:`the forward substitution <forward_substitution>`, I am left with the following system:

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

In the last row, I have

.. math::

   x_{n-1} = r_{n-1},

which has already been computed in :ref:`the forward substitution <forward_substitution>`.

Also, since I have

.. math::

   x_i = r_i - v_i x_{i+1},

I can compute :math:`x_i` one after another (sequentially from :math:`i = n-2` to :math:`i = 0`):

.. myliteralinclude:: /../../src/tdm.c
   :language: c
   :tag: backward substitution

Note again that ``q`` is shared among :math:`x_i` (output) and :math:`q_i` (input) in the code.

