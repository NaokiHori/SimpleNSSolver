
.. _linear_system:

######################
`src/linear_system.c`_
######################

.. _src/linear_system.c: https://github.com/NaokiHori/SimpleNSSolver/blob/main/src/linear_system.c

This file contains functions to solve the systems of linear equations:

.. math::

   A_{ij} x_j = b_i,

in particular it aims to solve this system when :math:`A_{ij}` is a tri-diagonal matrix.

As discussed in the temporal discretisations of :ref:`the momentum equation <smac_method>` and :ref:`the equation of internal energy <temperature_integration>`, when the diffusive terms are treated implicitly, I need to solve

.. math::

   u_i^*
   =
   u_i^k
   +
   \left( 1 - \frac{\gamma^k \Delta t}{2} \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{\delta^2}{\delta y^2} \right)^{-1}
   \left( 1 - \frac{\gamma^k \Delta t}{2} \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{\delta^2}{\delta x^2} \right)^{-1}
   \delta u_i

for the momentum equation, while

.. math::

   T^{k+1}
   =
   T^k
   +
   \left( 1 - \frac{\gamma^k \Delta t}{2} \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{\delta^2}{\delta y^2} \right)^{-1}
   \left( 1 - \frac{\gamma^k \Delta t}{2} \frac{\sqrt{Pr}}{\sqrt{Ra}} \frac{\delta^2}{\delta x^2} \right)^{-1}
   \delta T

for the temperature field, respectively.
Here, I only consider two-dimensional cases since the extension to three-dimensional cases is straightforward.

For the sake of notational simplicity, I write the right-hand-side of these equations as

.. math::

   \left( 1 - k \frac{\delta^2}{\delta y^2} \right)^{-1}
   \left( 1 - k \frac{\delta^2}{\delta x^2} \right)^{-1}
   q,

where :math:`k` is a constant value

.. math::

   \begin{cases}
      \text{Momentum}    & k \equiv \frac{\gamma^k \Delta t}{2} \frac{\sqrt{Pr}}{\sqrt{Ra}}, \\
      \text{Temperature} & k \equiv \frac{\gamma^k \Delta t}{2} \frac{1}{\sqrt{Pr} \sqrt{Ra}},
   \end{cases}

and :math:`q` is a scalar field which can be one of :math:`\ux`, :math:`\uy`, :math:`\uz` and :math:`T`.

.. mydetails:: Implementation (assigning to a variable ``prefactor``):

   .. myliteralinclude:: /../../src/fluid/integrate/ux.c
      :language: c
      :tag: gamma dt diffusivity / 2

   .. myliteralinclude:: /../../src/fluid/integrate/uy.c
      :language: c
      :tag: gamma dt diffusivity / 2

   .. myliteralinclude:: /../../src/fluid/integrate/uz.c
      :language: c
      :tag: gamma dt diffusivity / 2

   .. myliteralinclude:: /../../src/fluid/integrate/t.c
      :language: c
      :tag: gamma dt diffusivity / 2

   .. note::

      Although the grid sizes (:math:`\Delta x_{\xic}`, :math:`\Delta x_{\pic}`, :math:`\Delta y`, and :math:`\Delta z`) are constant in time, :math:`k` is time-dependent (:math:`\gamma^k \Delta t`) and thus I need to initialise the tri-diagonal coefficients every time step.

I look at the Laplace operator in the :math:`x` direction, i.e.

.. math::

   q^{\prime}
   =
   \left( 1 - k \frac{\delta^2}{\delta x^2} \right)^{-1}
   q,

or equivalently

.. math::

   \left( 1 - k \frac{\delta^2}{\delta x^2} \right)
   q^{\prime}
   =
   q.

Since I adopt the second-order-accurate central-difference scheme,

.. math::

   \frac{\delta^2}{\delta x^2}

is discretised as

.. math::

   \vat{\frac{\delta^2 q^{\prime}}{\delta x^2}}{\xic}
   \approx
     {\chi}_{\ximm} \vat{q^{\prime}}{\ximm}
   + {\chi}_{\xic } \vat{q^{\prime}}{\xic }
   + {\chi}_{\xipp} \vat{q^{\prime}}{\xipp},

where

.. math::

   \begin{aligned}
      {\chi}_{\ximm} & = \frac{1}{\Delta x_{\xic} \Delta x_{\xim}}, \\
      {\chi}_{\xipp} & = \frac{1}{\Delta x_{\xic} \Delta x_{\xip}}, \\
      {\chi}_{\xic } & = -{\chi}_{\ximm}-{\chi}_{\xipp}
   \end{aligned}

for :math:`q \leftarrow \ux`, otherwise

.. math::

   \vat{\frac{\delta^2 q^{\prime}}{\delta x^2}}{\pic}
   \approx
     {\chi}_{\pimm} \vat{q^{\prime}}{\pimm}
   + {\chi}_{\pic } \vat{q^{\prime}}{\pic }
   + {\chi}_{\pipp} \vat{q^{\prime}}{\pipp},

where

.. math::

   {\chi}_{\yimm} & = \frac{1}{\Delta x_{\yic} \Delta x_{\yim}}, \\
   {\chi}_{\yipp} & = \frac{1}{\Delta x_{\yic} \Delta x_{\yip}}, \\
   {\chi}_{\yic } & = -{\chi}_{\yimm}-{\chi}_{\yipp}.

Thus, I have

.. math::

   \left( 1 - k \frac{\delta^2}{\delta x^2} \right)
   \vat{q^{\prime}}{\xic}
   & \approx
              - k \vat{\chi}{\xipp} \vat{q^{\prime}}{\xipp}
   + \left( 1 - k \vat{\chi}{\xic } \vat{q^{\prime}}{\xic } \right)
              - k \vat{\chi}{\ximm} \vat{q^{\prime}}{\ximm} \\
   & \equiv
     \vat{u}{\xic} \vat{q^{\prime}}{\xipp}
   + \vat{c}{\xic} \vat{q^{\prime}}{\xic }
   + \vat{l}{\xic} \vat{q^{\prime}}{\ximm} \\
   & =
   \vat{q}{\xic}

for :math:`q \leftarrow \ux`, otherwise

.. math::

   \left( 1 - k \frac{\delta^2}{\delta x^2} \right)
   \vat{q^{\prime}}{\pic}
   & \approx
              - k \vat{\chi}{\pipp} \vat{q^{\prime}}{\pipp}
   + \left( 1 - k \vat{\chi}{\pic } \vat{q^{\prime}}{\pic } \right)
              - k \vat{\chi}{\pimm} \vat{q^{\prime}}{\pimm} \\
   & \equiv
     \vat{u}{\pic} \vat{q^{\prime}}{\pipp}
   + \vat{c}{\pic} \vat{q^{\prime}}{\pic }
   + \vat{l}{\pic} \vat{q^{\prime}}{\pimm} \\
   & =
   \vat{q}{\pic}.

I can write them using the following matrix forms (in particular the tri-diagonal linear systems):

.. math::

   \newcommand\ia{\frac{3}{2}}
   \newcommand\ib{\frac{5}{2}}
   \newcommand\id{i-\frac{1}{2}}
   \newcommand\ie{i+\frac{1}{2}}
   \newcommand\if{i+\frac{3}{2}}
   \newcommand\ih{\text{isize}-\frac{3}{2}}
   \newcommand\ii{\text{isize}-\frac{1}{2}}
   \begin{bmatrix}
      c_{\ia} & u_{\ia} & \cdots & 0       & 0       & 0       & \cdots & 0       & 0       \\
      l_{\ib} & c_{\ib} & \cdots & 0       & 0       & 0       & \cdots & 0       & 0       \\
      \vdots  & \vdots  & \ddots & \vdots  & \vdots  & \vdots  &        & \vdots  & \vdots  \\
      0       & 0       & \cdots & c_{\id} & u_{\id} & 0       & \cdots & 0       & 0       \\
      0       & 0       & \cdots & l_{\ie} & c_{\ie} & u_{\ie} & \cdots & 0       & 0       \\
      0       & 0       & \cdots & 0       & l_{\if} & c_{\if} & \cdots & 0       & 0       \\
      \vdots  & \vdots  &        & \vdots  & \vdots  & \vdots  & \ddots & \vdots  & \vdots  \\
      0       & 0       & \cdots & 0       & 0       & 0       & \cdots & c_{\ih} & u_{\ih} \\
      0       & 0       & \cdots & 0       & 0       & 0       & \cdots & l_{\ii} & c_{\ii}
   \end{bmatrix}
   \begin{bmatrix}
      \vat{q^{\prime}}{\ia} \\
      \vat{q^{\prime}}{\ib} \\
      \vdots                \\
      \vat{q^{\prime}}{\id} \\
      \vat{q^{\prime}}{\ie} \\
      \vat{q^{\prime}}{\if} \\
      \vdots                \\
      \vat{q^{\prime}}{\ih} \\
      \vat{q^{\prime}}{\ii}
   \end{bmatrix}
   =
   \begin{bmatrix}
      \vat{q}{\ia} \\
      \vat{q}{\ib} \\
      \vdots       \\
      \vat{q}{\id} \\
      \vat{q}{\ie} \\
      \vat{q}{\if} \\
      \vdots       \\
      \vat{q}{\ih} \\
      \vat{q}{\ii}
   \end{bmatrix}

for :math:`q \leftarrow \ux`, otherwise

.. math::

   \newcommand\ia{1}
   \newcommand\ib{2}
   \newcommand\id{i-1}
   \newcommand\ie{i  }
   \newcommand\if{i+1}
   \newcommand\ih{\text{isize}-1}
   \newcommand\ii{\text{isize}  }
   \begin{bmatrix}
      c_{\ia} & u_{\ia} & \cdots & 0       & 0       & 0       & \cdots & 0       & 0       \\
      l_{\ib} & c_{\ib} & \cdots & 0       & 0       & 0       & \cdots & 0       & 0       \\
      \vdots  & \vdots  & \ddots & \vdots  & \vdots  & \vdots  &        & \vdots  & \vdots  \\
      0       & 0       & \cdots & c_{\id} & u_{\id} & 0       & \cdots & 0       & 0       \\
      0       & 0       & \cdots & l_{\ie} & c_{\ie} & u_{\ie} & \cdots & 0       & 0       \\
      0       & 0       & \cdots & 0       & l_{\if} & c_{\if} & \cdots & 0       & 0       \\
      \vdots  & \vdots  &        & \vdots  & \vdots  & \vdots  & \ddots & \vdots  & \vdots  \\
      0       & 0       & \cdots & 0       & 0       & 0       & \cdots & c_{\ih} & u_{\ih} \\
      0       & 0       & \cdots & 0       & 0       & 0       & \cdots & l_{\ii} & c_{\ii}
   \end{bmatrix}
   \begin{bmatrix}
      \vat{q^{\prime}}{\ia} \\
      \vat{q^{\prime}}{\ib} \\
      \vdots                \\
      \vat{q^{\prime}}{\id} \\
      \vat{q^{\prime}}{\ie} \\
      \vat{q^{\prime}}{\if} \\
      \vdots                \\
      \vat{q^{\prime}}{\ih} \\
      \vat{q^{\prime}}{\ii}
   \end{bmatrix}
   =
   \begin{bmatrix}
      \vat{q}{\ia} \\
      \vat{q}{\ib} \\
      \vdots       \\
      \vat{q}{\id} \\
      \vat{q}{\ie} \\
      \vat{q}{\if} \\
      \vdots       \\
      \vat{q}{\ih} \\
      \vat{q}{\ii}
   \end{bmatrix}.

.. mydetails:: Boundary treatments

   I take :math:`\ux` as an example.
   Originally, in the vicinity of the walls (one-grid apart from the walls, at :math:`\frac{3}{2}` or at :math:`\text{isize}-\frac{1}{2}`), I have

   .. math::

      l_{\frac{3}{2}} q^{\prime}_{\text{left wall}}
      +
      c_{\frac{3}{2}} q^{\prime}_{\frac{3}{2}}
      +
      u_{\frac{3}{2}} q^{\prime}_{\frac{5}{2}}
      =
      q_{\frac{3}{2}},

   .. math::

      l_{\text{isize}-\frac{1}{2}} q^{\prime}_{\text{isize}-\frac{3}{2}}
      +
      c_{\text{isize}-\frac{1}{2}} q^{\prime}_{\text{isize}-\frac{1}{2}}
      +
      u_{\text{isize}-\frac{1}{2}} q^{\prime}_{\text{right wall}}
      =
      q_{\text{isize}-\frac{1}{2}}.

   Since :math:`q` is the delta form (:math:`q = \delta \ux \equiv \ux^* - \ux^k`), :math:`q` is equal to :math:`0` on the boundaries

   .. math::

      & q^{\prime}_{\text{left wall}}  \equiv 0, \\
      & q^{\prime}_{\text{right wall}} \equiv 0,

   because I assume that the walls are impermeable, and the Neumann boundary condition is imposed on the scalar potential :math:`\psi`, indicating that the correction step does not alter the values on the walls.
   Thus, I conclude

   .. math::

      c_{\frac{3}{2}} q^{\prime}_{\frac{3}{2}}
      +
      u_{\frac{3}{2}} q^{\prime}_{\frac{5}{2}}
      =
      q_{\frac{3}{2}},

   .. math::

      l_{\text{isize}-\frac{1}{2}} q^{\prime}_{\text{isize}-\frac{3}{2}}
      +
      c_{\text{isize}-\frac{1}{2}} q^{\prime}_{\text{isize}-\frac{1}{2}}
      =
      q_{\text{isize}-\frac{1}{2}},

   which are the corrections taking into account the boundary conditions.

Since :math:`y` and :math:`z` directions are homogeneous (periodic boundary conditions are imposed and grids are equidistantly placed), matrix in the left-hand side describing the Laplace operators in the :math:`y` and :math:`z` directions

.. math::

   \frac{\delta^2}{\delta y^2}, \frac{\delta^2}{\delta z^2}

are

.. math::

   \begin{bmatrix}
      c       & u       & \cdots & 0       & 0       & 0       & \cdots & 0       & l       \\
      l       & c       & \cdots & 0       & 0       & 0       & \cdots & 0       & 0       \\
      \vdots  & \vdots  & \ddots & \vdots  & \vdots  & \vdots  &        & \vdots  & \vdots  \\
      0       & 0       & \cdots & c       & u       & 0       & \cdots & 0       & 0       \\
      0       & 0       & \cdots & l       & c       & u       & \cdots & 0       & 0       \\
      0       & 0       & \cdots & 0       & l       & c       & \cdots & 0       & 0       \\
      \vdots  & \vdots  &        & \vdots  & \vdots  & \vdots  & \ddots & \vdots  & \vdots  \\
      0       & 0       & \cdots & 0       & 0       & 0       & \cdots & c       & u       \\
      u       & 0       & \cdots & 0       & 0       & 0       & \cdots & l       & c
   \end{bmatrix},

where the coefficients are

.. math::

   u = l &=   - k \frac{1}{\Delta y^2}, \\
   c     &= 1 + k \frac{2}{\Delta y^2},

or

.. math::

   u = l &=   - k \frac{1}{\Delta z^2}, \\
   c     &= 1 + k \frac{2}{\Delta z^2},

respectively.

The tri-diagonal matrix is solved by :ref:`the Thomas algorithm <tdm>`.

.. note::

   In this project, the solver requests that all systems are present on the memory, i.e. one process should know everything about the system to be solved.
   For systems in the :math:`x` direction, it is by default fulfilled since the domain is not decomposed.
   On the other hand, for systems in the :math:`y` and :math:`z` directions, this is not satisfied.
   Thus I need `the pencil rotations <https://github.com/NaokiHori/SimpleDecomp>`_ to solve linear systems in the :math:`y` and :math:`z` directions.

I need three steps to solve a linear system using this module:

 * Initialisation

 * Execution

 * Finalisation

In this file the initialisation and the finalisation are implemented.
The execution step is not handled since it is just to call :ref:`the Thomas algorithm <tdm>`.

**************
Initialisation
**************

.. mydeclare:: /../../src/linear_system.c
   :language: c
   :tag: linear_system_init

This function plays the following roles:

   #. allocating the internal buffers (if needed),

   #. preparing for the pencil rotations (if needed),

   #. initialising (the coefficients of) the tri-diagonal matrices (if needed).

================
Allocate buffers
================

``y1pencil`` and ``z2pencil`` are only needed when the corresponding direction is treated implicitly.

.. myliteralinclude:: /../../src/linear_system.c
   :language: c
   :tag: allocate pencils if needed

===========================
Initialise pencil rotations
===========================

``xy`` rotations and ``xz`` rotations are only needed when the :math:`y` and the :math:`z` directions are treated implicitly, respectively.

.. myliteralinclude:: /../../src/linear_system.c
   :language: c
   :tag: between x1 and y1

.. myliteralinclude:: /../../src/linear_system.c
   :language: c
   :tag: between x1 and z2

=======================
Initialise the matrices
=======================

Solver in each direction is only needed when the direction is treated implicitly.

.. myliteralinclude:: /../../src/linear_system.c
   :language: c
   :tag: Thomas algorithm in x direction

.. myliteralinclude:: /../../src/linear_system.c
   :language: c
   :tag: Thomas algorithm in y direction

.. myliteralinclude:: /../../src/linear_system.c
   :language: c
   :tag: Thomas algorithm in z direction

