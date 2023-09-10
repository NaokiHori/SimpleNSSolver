
.. _incompressibility_and_poisson_equation:

######################################
Incompressibility and Poisson equation
######################################

********
Overview
********

In this part, I consider to discretise the incompressibility constraint

.. math::

   \der{u_i}{x_i} = 0

at the cell center :math:`\left( \pic, \pjc \right)`.

Also the discretisation of the Poisson equation

.. math::

   \der{}{x_i} \der{\psi}{x_i}
   =
   \frac{1}{\Delta t} \der{u_i^*}{x_i}

is derived.

.. seealso::

   Please refer to :ref:`the SMAC method <smac_method>` to see why the Poisson equation appears.

*************************************
Discrete incompressibility constraint
*************************************

The incompressibility constraint at cell center :math:`\left( \pic, \pjc \right)` is discretised as

.. math::

   \dder{\ux}{x}
   +
   \dder{\uy}{y}
   =
   0.

.. mydetails:: Derivations

   As discussed in :ref:`the grid arrangement <domain_setup>`, in this project, cell centers (where the pressure and the temperature are defined) are positioned in the middle of the surrounding cell faces (where the velocities are defined), i.e.

   .. math::

      \vat{x}{\pic}
      & =
      \frac{1}{2} \vat{x}{\pip}
      +
      \frac{1}{2} \vat{x}{\pim}, \\
      \vat{y}{\pjc}
      & =
      \frac{1}{2} \vat{y}{\pjp}
      +
      \frac{1}{2} \vat{y}{\pjm}.

   Since the velocities are defined at the cell faces, it is natural to describe the incompressibility constraint at cell centers:

   .. math::

      \der{u_i}{x_i} \left( x, y \right)
      \approx
      \vat{\dder{u_i}{x_i}}{\pic,\pjc}
      =
      \frac{\vat{\ux}{\pip} - \vat{\ux}{\pim}}{\vat{x}{\pip} - \vat{x}{\pim}}
      +
      \frac{\vat{\uy}{\pjp} - \vat{\uy}{\pjm}}{\vat{y}{\pjp} - \vat{y}{\pjm}}
      =
      0,

   giving the discrete incompressibility constraint defined at the cell centers.

   .. note::

      The continuity is *only* valid at each cell center.
      This is clearly different from the continuous domain, where all infinitesimal control volumes obey it.

*************************
Discrete Poisson equation
*************************

The discrete Poisson equation at :math:`\left( \pic, \pjc \right)` is given by

.. math::

   \dder{}{x} \left(
      \dder{\psi}{x}
   \right)
   +
   \dder{}{y} \left(
      \dder{\psi}{y}
   \right)
   =
   \cdots.

.. mydetails:: Derivations

   As discussed in :ref:`the temporal discretisation <temporal_discretisation>`, the continuity is enforced by the correction:

   .. math::

      u_i^{n+1}
      =
      u_i^{*}
      -
      \der{\psi}{x_i} \Delta t.

   This relation holds on the cell faces:

   .. math::

      \vat{\ux^{n+1}}{\pip} &= \vat{\ux^*}{\pip} - \vat{\der{\psi}{x}}{\pip} \Delta t, \\
      \vat{\ux^{n+1}}{\pim} &= \vat{\ux^*}{\pim} - \vat{\der{\psi}{x}}{\pim} \Delta t, \\
      \vat{\uy^{n+1}}{\pjp} &= \vat{\uy^*}{\pjp} - \vat{\der{\psi}{y}}{\pjp} \Delta t, \\
      \vat{\uy^{n+1}}{\pjm} &= \vat{\uy^*}{\pjm} - \vat{\der{\psi}{y}}{\pjm} \Delta t,

   where :math:`\psi` is a scalar potential which projects the non-solenoidal velocity field to the solenoidal one.
   Assigning this four equations to the discrete incompressibility constraint, I have

   .. math::

      \frac{
          \vat{\der{\psi}{x}}{\pip}
        - \vat{\der{\psi}{x}}{\pim}
      }{
          \vat{x}{\pip}
        - \vat{x}{\pim}
      }
      +
      \frac{
          \vat{\der{\psi}{y}}{\pjp}
        - \vat{\der{\psi}{y}}{\pjm}
      }{
          \vat{y}{\pjp}
        - \vat{y}{\pjm}
      }
      =
      \frac{1}{\Delta t}
      \frac{\vat{\ux^*}{\pip} - \vat{\ux^*}{\pim}}{\vat{x}{\pip} - \vat{x}{\pim}}
      +
      \frac{\vat{\uy^*}{\pjp} - \vat{\uy^*}{\pjm}}{\vat{y}{\pjp} - \vat{y}{\pjm}},

   which is the discrete pressure equation

   .. math::

      \dder{}{x_i} \dder{\psi}{x_i}
      =
      \frac{1}{\Delta t} \dder{u_i^*}{x_i}.

   .. note::

      When the grid sizes are not uniform, the above discretisation is different from the second-order-accurate approximation of the Taylor series expansion of

      .. math::

         \der{}{x_i} \der{\psi}{x_i}.

      If one uses this, the incompressibility constraint is not satisfied.

