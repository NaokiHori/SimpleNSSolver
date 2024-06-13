
.. _momentum_integration:

.. include:: /references.txt

################
Momentum Balance
################

The momentum balance is

.. math::

    \pder{u_i}{t}
    =
    P_i
    +
    A_i
    +
    D_{i1}
    +
    D_{i2}
    +
    D_{i3},

where :math:`P_i` and :math:`A_i` are the pressure-gradient terms:

.. math::

    P_i
    \equiv
    -
    \dmompre{i},

and the advective terms:

.. math::

    A_i
    \equiv
    -
    \dmomadv{i}{1}
    -
    \dmomadv{i}{2}
    -
    \dmomadv{i}{3},

respectively.

:math:`D_{ij}` is the diffusive term involving spatial differentiation in the :math:`j`-th direction:

.. math::

    D_{ij}
    \equiv
    \dmomdif{j}{i}
    \,\,
    (\text{no summation over}\,j).

See :ref:`the spatial discretization <discrete_momentum_balance>`.

The temporal discretization for each Runge-Kutta iteration leads to

.. math::

    \Delta u_i
    &
    =
    \gamma^k \Delta t P_i^k
    +
    \alpha^k \Delta t \left( A_i^{k  } + D_i^{k  } \right)
    +
    \beta^k  \Delta t \left( A_i^{k-1} + D_i^{k-1} \right),

    u_i^*
    &
    =
    u_i^k
    +
    \Delta u_i,

when all advective and diffusive terms are treated explicitly, while

.. math::

    \newcommand{\lap}[2]{
        {#2} \gamma^k \Delta t
        \frac{\sqrt{Pr}}{\sqrt{Ra}}
        \frac{1}{J}
        \dif{}{\gcs{#1}}
        \left(
            \frac{J}{\sfact{#1}}
            \frac{1}{\sfact{#1}}
            \dif{}{\gcs{#1}}
        \right)
    }
    \Delta u_i
    &
    =
    \gamma^k \Delta t P_i^k
    +
    \alpha^k \Delta t A_i^{k  }
    +
    \beta^k  \Delta t A_i^{k-1}
    +
    \gamma^k \Delta t \left( D_{i1}^k + D_{i2}^k + D_{i3}^k \right),

    u_i^*
    &
    =
    u_i^k
    +
    \left\{
        1
        -
        \lap{3}{c}
    \right\}^{-1}
    \left\{
        1
        -
        \lap{2}{c}
    \right\}^{-1}
    \left\{
        1
        -
        \lap{1}{c}
    \right\}^{-1}
    \Delta u_i,

when diffusive terms are treated implicitly.

Although we obtain a new velocity field, this does not satisfy :ref:`the incompressibility <discrete_incompressibility>` in general, which necessitates the additional procedure below.

***********
SMAC method
***********

In addition to the incompressibility constraint, we need to somehow update the pressure field as well, which we do not have any equation such as:

.. math::

    \pder{p}{t}
    =
    \cdots.

To overcome these issues, we adopt the Simplified Marker And Cell (SMAC) method (|AMSDEN1970|), which is a two-step method.

===============
Prediction Step
===============

In the first step (prediction step), momentum equation is integrated in time, without taking care of the mass conservation, as discussed above.
Basically the procedure is identical to :ref:`how we handle the temperature field <temperature_integration>`.
First, the explicit and implicit terms are calculated and stored to the corresponding buffers:

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: compute right-hand-side terms, which are added to buffers

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: compute right-hand-side terms, which are added to buffers

.. myliteralinclude:: /../../src/fluid/predict/uz.c
    :language: c
    :tag: compute right-hand-side terms, which are added to buffers

The stored values are used to compute :math:`\Delta u_i`:

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: compute increments

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: compute increments

.. myliteralinclude:: /../../src/fluid/predict/uz.c
    :language: c
    :tag: compute increments

When necessary, linear systems are solved to take care of the implicit treatments:

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: solve linear systems in x

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: solve linear systems in y

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: solve linear systems in z

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: solve linear systems in x

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: solve linear systems in y

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: solve linear systems in z

.. myliteralinclude:: /../../src/fluid/predict/uz.c
    :language: c
    :tag: solve linear systems in x

.. myliteralinclude:: /../../src/fluid/predict/uz.c
    :language: c
    :tag: solve linear systems in y

.. myliteralinclude:: /../../src/fluid/predict/uz.c
    :language: c
    :tag: solve linear systems in z

Finally the velocity field is updated:

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: update velocity field

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: update velocity field

.. myliteralinclude:: /../../src/fluid/predict/uz.c
    :language: c
    :tag: update velocity field

===============
Correction Step
===============

The updated velocity field :math:`u_i^*`, which in general violates the incompressibility, is corrected in the second step (correction step).
The idea is mathematically written as

.. math::

    u_i^{k+1}
    =
    u_i^*
    -
    \gamma^k \Delta t
    \frac{1}{\sfact{i}}
    \dif{\psi}{\gcs{i}},

where :math:`\psi` is a scalar potential to be given.
By taking the (discrete) divergence, we obtain

.. math::

    \frac{1}{J}
    \dif{}{\gcs{i}}
    \left(
        \frac{J}{\sfact{i}}
        \vel{i}^{k+1}
    \right)
    =
    \frac{1}{J}
    \dif{}{\gcs{i}}
    \left(
        \frac{J}{\sfact{i}}
        \vel{i}^*
    \right)
    -
    \gamma^k \Delta t
    \frac{1}{J}
    \dif{}{\gcs{i}}
    \left(
        \frac{J}{\sfact{i}}
        \frac{1}{\sfact{i}}
        \dif{\psi}{\gcs{i}}
    \right).

By requesting the (discrete) incompressibility constraint on the new velocity field (namely the left-hand-side term to be zero), we obtain

.. math::

    \frac{1}{J}
    \dif{}{\gcs{i}}
    \left(
        \frac{J}{\sfact{i}}
        \frac{1}{\sfact{i}}
        \dif{\psi}{\gcs{i}}
    \right)
    =
    \frac{1}{\gamma^k \Delta t}
    \frac{1}{J}
    \dif{}{\gcs{i}}
    \left(
        \frac{J}{\sfact{i}}
        \vel{i}^*
    \right),

which is a Poisson equation with respect to :math:`\psi`.

The right-hand-side term is computed as follows in the code:

.. myliteralinclude:: /../../src/fluid/compute_potential.c
    :language: c
    :tag: compute right-hand side of Poisson equation

After :ref:`solving the Poisson equation <poisson_equation>`, the velocity field is corrected as follows in the code:

.. myliteralinclude:: /../../src/fluid/correct/ux.c
    :language: c
    :tag: correct x velocity

.. myliteralinclude:: /../../src/fluid/correct/uy.c
    :language: c
    :tag: correct y velocity

.. myliteralinclude:: /../../src/fluid/correct/uz.c
    :language: c
    :tag: correct z velocity

=============================
Pressure and Scalar Potential
=============================

Finally we relate the pressure field with the scalar potential to close the system.
Each Runge-Kutta step when all diffusive terms are treated explicitly is given by

.. math::

    \Delta u_i
    &
    =
    \gamma^k \Delta t P_i^k
    +
    \alpha^k \Delta t \left( A_i^{k  } + D_{i1}^{k  } + D_{i2}^{k  } + D_{i3}^{k  } \right)
    +
    \beta^k  \Delta t \left( A_i^{k-1} + D_{i1}^{k-1} + D_{i2}^{k-1} + D_{i3}^{k-1} \right),

    u_i^*
    &
    =
    u_i^k
    +
    \Delta u_i,

    u_i^{k+1}
    &
    =
    u_i^*
    -
    \gamma^k \Delta t
    \frac{1}{\sfact{i}}
    \dif{\psi}{\gcs{i}}.

Summing all three steps yield

.. math::

    u_i^{k+1}
    =
    u_i^k
    -
    \gamma^k \Delta t
    \frac{1}{\sfact{i}}
    \dif{\psi}{\gcs{i}}
    +
    \gamma^k \Delta t P_i^k
    +
    \alpha^k \Delta t \left( A_i^{k  } + D_{i1}^{k  } + D_{i2}^{k  } + D_{i3}^{k  } \right)
    +
    \beta^k  \Delta t \left( A_i^{k-1} + D_{i1}^{k-1} + D_{i2}^{k-1} + D_{i3}^{k-1} \right).

Since :ref:`the pressure field should be treated implicitly in time <implicit_treatment>`, we find a requirement:

.. math::

    -
    \gamma^k \Delta t
    \frac{1}{\sfact{i}}
    \dif{\psi}{\gcs{i}}
    +
    \gamma^k \Delta t P_i^k
    =
    \gamma^k \Delta t P_i^{k+1},

or equivalently

.. math::

    p^{k+1}
    =
    p^k
    +
    \psi.

Next, we consider cases where the diffusive terms are partially treated implicitly in time; for instance, the following relation holds when the diffusive terms are treated implicitly only in :math:`x` direction:

.. math::

    \Delta u_i
    &
    =
    \gamma^k \Delta t P_i^k
    +
    \alpha^k \Delta t \left( A_i^{k  } + D_{i2}^{k  } + D_{i3}^{k  } \right)
    +
    \beta^k  \Delta t \left( A_i^{k-1} + D_{i2}^{k-1} + D_{i3}^{k-1} \right)
    +
    \gamma^k \Delta t D_{i1}^{k  },

    u_i^*
    &
    =
    u_i^k
    +
    \left\{
        1
        -
        \lap{1}{c}
    \right\}^{-1}
    \Delta u_i,

    u_i^{k+1}
    &
    =
    u_i^*
    -
    \gamma^k \Delta t
    \frac{1}{\sfact{i}}
    \dif{\psi}{\gcs{i}}.

Note that, in the second step, the left-hand side is :math:`u_i^*` while it should be :math:`u_i^{k+1}` but is unknown.
By requesting the pressure-gradient term to be implicit in time and with some algebra, we obtain

.. math::

    \gamma^k \Delta t
    \left\{
        1
        -
        \lap{1}{c}
    \right\}
    \frac{1}{\sfact{i}}
    \dif{\psi}{\gcs{i}}
    -
    \gamma^k \Delta t
    P_i^k
    =
    -
    \gamma^k \Delta t
    P_i^{k+1},

or equivalently

.. math::

    \frac{1}{\sfact{i}}
    \dif{p^{k+1}}{\gcs{i}}
    =
    \frac{1}{\sfact{i}}
    \dif{p^k}{\gcs{i}}
    +
    \left\{
        1
        -
        \lap{1}{c}
    \right\}
    \frac{1}{\sfact{i}}
    \dif{\psi}{\gcs{i}}.

Since the spatial-differential operators are interchangeable, we simplify the relation to

.. math::

    p^{k+1}
    =
    p^k
    +
    \psi
    -
    \lap{1}{c}
    \psi.

Similarly, when all diffusive terms are treated implicitly, we have

.. math::

    \Delta u_i
    &
    =
    \gamma^k \Delta t P_i^k
    +
    \alpha^k \Delta t A_i^{k  }
    +
    \beta^k  \Delta t A_i^{k-1}
    +
    \gamma^k \Delta t \left( D_{i1}^{k  } + D_{i2}^{k  } + D_{i3}^{k  } \right),

    u_i^*
    &
    =
    u_i^k
    +
    \left\{
        1
        -
        \lap{3}{c}
    \right\}^{-1}
    \left\{
        1
        -
        \lap{2}{c}
    \right\}^{-1}
    \left\{
        1
        -
        \lap{1}{c}
    \right\}^{-1}
    \Delta u_i,

    u_i^{k+1}
    &
    =
    u_i^*
    -
    \gamma^k \Delta t
    \frac{1}{\sfact{i}}
    \dif{\psi}{\gcs{i}},

    p^{k+1}
    &
    =
    p^k
    +
    \psi
    -
    \lap{1}{c}
    \psi
    -
    \lap{2}{c}
    \psi
    -
    \lap{3}{c}
    \psi.

Updating pressure field using the scalar potential is implemented as follows:

.. myliteralinclude:: /../../src/fluid/update_pressure.c
    :language: c
    :tag: update pressure field using scalar potential

The explicit contribution, which is always present, is given here:

.. myliteralinclude:: /../../src/fluid/update_pressure.c
    :language: c
    :tag: explicit contribution

The implicit contributions, which is needed when the Laplace operator in the direction is implicitly treated, are given here:

.. myliteralinclude:: /../../src/fluid/update_pressure.c
    :language: c
    :tag: x implicit contribution

.. myliteralinclude:: /../../src/fluid/update_pressure.c
    :language: c
    :tag: y implicit contribution

.. myliteralinclude:: /../../src/fluid/update_pressure.c
    :language: c
    :tag: z implicit contribution

