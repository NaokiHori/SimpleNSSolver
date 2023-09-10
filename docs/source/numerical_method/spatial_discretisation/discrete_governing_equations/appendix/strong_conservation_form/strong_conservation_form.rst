########################
Strong conservation form
########################

I consider to project the governing equations to a new coordinate system.

To begin with, I take the incompressibility constraint:

.. math::

   \der{u_i}{x_i}
   =
   0

as an example.
I try to replace the derivatives :math:`\partial / \partial x_i` by :math:`\partial / \partial \xi^i` as

.. math::

   \der{\xi^j}{x_i} \der{u_i}{\xi^j}
   = \der{ \xi}{x} \der{\ux}{ \xi}
   + \der{\eta}{x} \der{\ux}{\eta}
   + \der{ \xi}{y} \der{\uy}{ \xi}
   + \der{\eta}{y} \der{\uy}{\eta}.

I would like to make the right-hand-side terms conservative, i.e. :math:`\partial / \partial \xi^j \left( \cdots \right)`.
To do so, I consider the relation between the old and the new coordinate systems:

.. math::

   \begin{pmatrix}
      \der{}{x} \\
      \der{}{y}
   \end{pmatrix}
   =
   \begin{pmatrix}
      \der{\xi}{x} & \der{\eta}{x} \\
      \der{\xi}{y} & \der{\eta}{y}
   \end{pmatrix}
   \begin{pmatrix}
      \der{}{ \xi} \\
      \der{}{\eta}
   \end{pmatrix}

or

.. math::

   \begin{pmatrix}
      \der{}{ \xi} \\
      \der{}{\eta}
   \end{pmatrix}
   =
   \begin{pmatrix}
      \der{x}{ \xi} & \der{y}{ \xi} \\
      \der{x}{\eta} & \der{y}{\eta}
   \end{pmatrix}
   \begin{pmatrix}
      \der{}{x} \\
      \der{}{y}
   \end{pmatrix}.

Since the second equation yields

.. math::

   \begin{pmatrix}
      \der{}{x} \\
      \der{}{y}
   \end{pmatrix}
   =
   \frac{1}{\der{x}{ \xi} \der{y}{\eta} - \der{y}{ \xi} \der{x}{\eta}}
   \begin{pmatrix}
      \der{y}{\eta} & -\der{x}{\eta} \\
      -
      \der{y}{ \xi} &  \der{x}{ \xi}
   \end{pmatrix}
   \begin{pmatrix}
      \der{}{ \xi} \\
      \der{}{\eta}
   \end{pmatrix},

I find

.. math::

   \frac{1}{
      \der{x}{ \xi} \der{y}{\eta}
      -
      \der{x}{\eta} \der{y}{ \xi}
   }
   \begin{pmatrix}
      \der{y}{\eta} & -\der{x}{\eta} \\
     -\der{y}{ \xi} &  \der{x}{ \xi}
   \end{pmatrix}
   =
   \begin{pmatrix}
     \der{ \xi}{x} & \der{ \xi}{y} \\
     \der{\eta}{x} & \der{\eta}{y}
   \end{pmatrix},

where the denominator is the determinant of the Jacobian matrix.

As a result, the incompressibility constraint leads to

.. math::

   \der{\xi^j}{x_i} \der{u_i}{\xi^j}
   &= \der{   x}{ \xi} \der{\ux}{ \xi}
    + \der{   x}{\eta} \der{\ux}{\eta}
    + \der{   y}{ \xi} \der{\uy}{ \xi}
    + \der{   y}{\eta} \der{\uy}{\eta} \\
   &= \frac{1}{J} \der{   y}{\eta} \der{\ux}{ \xi}
    - \frac{1}{J} \der{   y}{ \xi} \der{\ux}{\eta}
    - \frac{1}{J} \der{   x}{\eta} \der{\uy}{ \xi}
    + \frac{1}{J} \der{   x}{ \xi} \der{\uy}{\eta}.

Next, I consider to update the velocity, i.e. trying to write down the equation using the contra-variant component :math:`U^i`.
Assigning

.. math::

   \begin{pmatrix}
     \ux \\
     \uy
   \end{pmatrix}
   \equiv
   \begin{pmatrix}
     \frac{\partial x}{\partial t} \\
     \frac{\partial y}{\partial t}
   \end{pmatrix}
   =
   \begin{pmatrix}
     \der{   x}{ \xi} & \der{   x}{\eta} \\
     \der{   y}{ \xi} & \der{   y}{\eta}
   \end{pmatrix}
   \begin{pmatrix}
     \frac{\partial  \xi}{\partial t} \\
     \frac{\partial \eta}{\partial t}
   \end{pmatrix}
   =
   \begin{pmatrix}
     \der{   x}{ \xi} & \der{   x}{\eta} \\
     \der{   y}{ \xi} & \der{   y}{\eta}
   \end{pmatrix}
   \begin{pmatrix}
     U^x \\
     U^y
   \end{pmatrix}

yields

.. math::

   &+ \frac{1}{J} \der{   y}{\eta} \frac{\partial}{\partial  \xi} \left( \der{   x}{ \xi} U^x + \der{   x}{\eta} U^y \right) \\
   &- \frac{1}{J} \der{   y}{ \xi} \frac{\partial}{\partial \eta} \left( \der{   x}{ \xi} U^x + \der{   x}{\eta} U^y \right) \\
   &- \frac{1}{J} \der{   x}{\eta} \frac{\partial}{\partial  \xi} \left( \der{   y}{ \xi} U^x + \der{   y}{\eta} U^y \right) \\
   &+ \frac{1}{J} \der{   x}{ \xi} \frac{\partial}{\partial \eta} \left( \der{   y}{ \xi} U^x + \der{   y}{\eta} U^y \right).

Terms including :math:`U` are

.. math::

   \begin{aligned}
     &+ \frac{1}{J} \der{   y}{\eta} \frac{\partial}{\partial  \xi} \left( \der{   x}{ \xi} U^x \right) \\
     &- \frac{1}{J} \der{   y}{ \xi} \frac{\partial}{\partial \eta} \left( \der{   x}{ \xi} U^x \right) \\
     &- \frac{1}{J} \der{   x}{\eta} \frac{\partial}{\partial  \xi} \left( \der{   y}{ \xi} U^x \right) \\
     &+ \frac{1}{J} \der{   x}{ \xi} \frac{\partial}{\partial \eta} \left( \der{   y}{ \xi} U^x \right)
   \end{aligned}
   =
   \begin{aligned}
     &+ \frac{U^x}{J} \der{   y}{\eta} \frac{\partial}{\partial  \xi} \left( \der{   x}{ \xi} \right) + \frac{1}{J} \der{   y}{\eta} \frac{\partial U^x}{\partial  \xi} \der{   x}{ \xi} \\
     &- \frac{U^x}{J} \der{   y}{ \xi} \frac{\partial}{\partial \eta} \left( \der{   x}{ \xi} \right) - \frac{1}{J} \der{   y}{ \xi} \frac{\partial U^x}{\partial \eta} \der{   x}{ \xi} \\
     &- \frac{U^x}{J} \der{   x}{\eta} \frac{\partial}{\partial  \xi} \left( \der{   y}{ \xi} \right) - \frac{1}{J} \der{   x}{\eta} \frac{\partial U^x}{\partial  \xi} \der{   y}{ \xi} \\
     &+ \frac{U^x}{J} \der{   x}{ \xi} \frac{\partial}{\partial \eta} \left( \der{   y}{ \xi} \right) + \frac{1}{J} \der{   x}{ \xi} \frac{\partial U^x}{\partial \eta} \der{   y}{ \xi}
   \end{aligned}
   \equiv
   \begin{aligned}
     &\left( 1-1 \right) + \left( 1-2 \right) \\
     &\left( 2-1 \right) + \left( 2-2 \right) \\
     &\left( 3-1 \right) + \left( 3-2 \right) \\
     &\left( 4-1 \right) + \left( 4-2 \right)
   \end{aligned},

.. math::

   \left( 1-2 \right) + \left( 3-2 \right) =
   + \frac{1}{J} \der{   y}{\eta} \frac{\partial U^x}{\partial  \xi} \der{   x}{ \xi}
   - \frac{1}{J} \der{   x}{\eta} \frac{\partial U^x}{\partial  \xi} \der{   y}{ \xi}
   = \frac{1}{J} \left( \der{   x}{ \xi} \der{   y}{\eta} - \der{   x}{\eta} \der{   y}{ \xi} \right) \frac{\partial U^x}{\partial \xi}
   = \frac{1}{J} J \frac{\partial U^x}{\partial \xi},

.. math::

   \left( 2-2 \right) + \left( 4-2 \right) =
   - \frac{1}{J} \der{   y}{ \xi} \frac{\partial U^x}{\partial \eta} \der{   x}{ \xi}
   + \frac{1}{J} \der{   x}{ \xi} \frac{\partial U^x}{\partial \eta} \der{   y}{ \xi}
   = 0,

.. math::

   \begin{aligned}
     & \left( 1-1 \right) + \left( 2-1 \right) + \left( 3-1 \right) + \left( 4-1 \right)  \\
     & =
     + \frac{U^x}{J} \der{   y}{\eta} \frac{\partial}{\partial  \xi} \left( \der{   x}{ \xi} \right)
     - \frac{U^x}{J} \der{   y}{ \xi} \frac{\partial}{\partial \eta} \left( \der{   x}{ \xi} \right)
     - \frac{U^x}{J} \der{   x}{\eta} \frac{\partial}{\partial  \xi} \left( \der{   y}{ \xi} \right)
     + \frac{U^x}{J} \der{   x}{ \xi} \frac{\partial}{\partial \eta} \left( \der{   y}{ \xi} \right) \\
     & =
     + \frac{U^x}{J} \der{   y}{\eta} \frac{\partial}{\partial  \xi} \left( \der{   x}{ \xi} \right)
     - \frac{U^x}{J} \der{   y}{ \xi} \frac{\partial}{\partial  \xi} \left( \der{   x}{\eta} \right)
     - \frac{U^x}{J} \der{   x}{\eta} \frac{\partial}{\partial  \xi} \left( \der{   y}{ \xi} \right)
     + \frac{U^x}{J} \der{   x}{ \xi} \frac{\partial}{\partial  \xi} \left( \der{   y}{\eta} \right) \\
     & =
     \frac{U^x}{J} \left[
       + \der{   x}{ \xi} \frac{\partial}{\partial  \xi} \left( \der{   y}{\eta} \right)
       + \der{   y}{\eta} \frac{\partial}{\partial  \xi} \left( \der{   x}{ \xi} \right)
       - \der{   x}{\eta} \frac{\partial}{\partial  \xi} \left( \der{   y}{ \xi} \right)
       - \der{   y}{ \xi} \frac{\partial}{\partial  \xi} \left( \der{   x}{\eta} \right)
     \right] \\
     & =
     \frac{U^x}{J} \left[
       + \frac{\partial}{\partial  \xi} \left(
         \der{   x}{ \xi} \der{   y}{\eta}
       \right)
       - \frac{\partial}{\partial  \xi} \left(
         \der{   x}{\eta} \der{   y}{ \xi}
       \right)
     \right] \\
     & =
     \frac{U^x}{J} \frac{\partial}{\partial \xi} \left( \der{   x}{ \xi} \der{   y}{\eta} - \der{   x}{\eta} \der{   y}{ \xi} \right)
     = \frac{U^x}{J} \frac{\partial J}{\partial \xi},
   \end{aligned}

and as a result

.. math::

   \frac{1}{J} \frac{\partial}{\partial \xi} \left( J U^x \right).

Similarly I notice the terms involving :math:`U^y` yield

.. math::

   \frac{1}{J} \frac{\partial}{\partial \eta} \left( J U^y \right).

Thus I obtain

.. math::

   \frac{\partial u_i}{\partial x_i}
   =
   \frac{\partial \xi^j}{\partial x_i} \frac{\partial u_i}{\partial \xi^j}
   =
   \frac{1}{J} \frac{\partial}{\partial \xi^i} \left( J U^i \right)
   =
   0.

Since

.. math::

   \begin{pmatrix}
     U^x \\
     U^y
   \end{pmatrix}
   =
   \frac{1}{J}
   \begin{pmatrix}
      \der{   y}{\eta} & -\der{   x}{\eta} \\
     -\der{   y}{ \xi} &  \der{   x}{ \xi}
   \end{pmatrix}
   \begin{pmatrix}
     \ux \\
     \uy
   \end{pmatrix}
   =
   \begin{pmatrix}
     \der{   x}{ \xi} & \der{   y}{ \xi} \\
     \der{   x}{\eta} & \der{   y}{\eta}
   \end{pmatrix}
   \begin{pmatrix}
     \ux \\
     \uy
   \end{pmatrix}

or

.. math::

   U^j = \frac{\partial \xi^j}{\partial x_i} u_i,

I finally obtain

.. math::

   \frac{\partial u_i}{\partial x_i}
   =
   \frac{\partial \xi^j}{\partial x_i} \frac{\partial u_i}{\partial \xi^j}
   =
   \frac{1}{J} \frac{\partial}{\partial \xi^j} \left( J \frac{\partial \xi^j}{\partial x_i} u_i \right)
   =
   0,

which is the conservative form of the incompressibility constraint defined in the new coordinate system (strong conservation form).

By comparing the second term with the third one, I notice

.. math::

   \frac{\partial \xi^j}{\partial x_i} \frac{\partial u_i}{\partial \xi^j}
   =
   \frac{1}{J} \left\{ \frac{\partial}{\partial \xi^j} \left( J \frac{\partial \xi^j}{\partial x_i} \right) \right\} u_i
   + \frac{1}{J} J \frac{\partial \xi^j}{\partial x_i} \frac{\partial u_i}{\partial \xi^j}
   =
   \frac{1}{J} \left\{ \frac{\partial}{\partial \xi^j} \left( J \frac{\partial \xi^j}{\partial x_i} \right) \right\} u_i
   + \frac{\partial \xi^j}{\partial x_i} \frac{\partial u_i}{\partial \xi^j}
   =
   0,

and thus I find the following important identity:

.. math::

   \frac{1}{J} \left\{ \frac{\partial}{\partial \xi^j} \left( J \frac{\partial \xi^j}{\partial x_i} \right) \right\} u_i = 0.

By following the same procedure, I notice that this relation even holds for general tensors (not limited to the first-order tensors).
This relation is used to derive the strong conservation forms of the momentum and internal energy balances.

.. note::

   Since the current projection keeps the orthogonality, the Jacobian matrix is a diagonal matrix.
   Thus, :math:`\der{\xi^j}{x_i}` can be reduced to :math:`\der{\xi^i}{x_i}`.

