
.. _implicit_treatment:

.. include:: /references.txt

##################
Implicit Treatment
##################

**************************************
Time Scales and Time-Step Restrictions
**************************************

In the Navier-Stokes equations, there are advective, pressure-gradient, and diffusive terms, each with different time scales.
Numerically, when a term is treated explicitly in time, information associated with it should not propagate a distance greater than the grid size in one time step.
In other words, to allow the information to travel further, the term should be treated implicitly.
Here, we elaborate on the three terms and their time scales to understand their effects on the overall time-marching process.

#. Advective terms

    By using the fluid velocity and the grid size as the reference velocity and length scales, we find that the advective terms impose a constraint:

    .. math::

        &
        \Delta t_{adv}
        =
        C
        \frac{\sfact{i}}{\vel{i}},

        &
        C
        <
        1,

    where :math:`C` is a non-dimensional number known as the Courant number.

    Recall that the grid sizes are unity in the computational coordinate system, while the velocity is divided by the scale factor.
    See :ref:`the equations in strong conservation forms <strong_conservation_form>`.

    Due to advective effects, :math:`\Delta t` should be reduced as the resolution becomes finer, being proportional to the spatial resolution.
    Since the advective terms are non-linear, treating them implicitly is not straightforward, thus this constraint is always present in this project.

#. Diffusive terms

    By adopting the diffusivities and the grid size as reference scales, we find that the diffusive terms impose another constraint:

    .. math::

        &
        \Delta t_{dif}
        =
        F
        \min
        \left(
            \frac{\sqrt{Ra}}{\sqrt{Pr}} ,
            \sqrt{Pr} \sqrt{Ra}
        \right)
        \sfact{i}^2,

        &
        F
        <
        1,

    where :math:`F` is a non-dimensional number known as the Fourier number.
    Note that the momentum and temperature fields have different diffusivities.

    As the spatial resolution is refined, :math:`\Delta t` should be reduced quadratically.
    This criterion can make computational costs prohibitive, especially for wall-bounded turbulent flows where the wall-normal grid sizes must be extremely small close to the walls to resolve boundary layers.

    Since the diffusive terms are linear, this restriction can be eliminated by treating them implicitly, as elaborated later on this page.

#. Pressure-gradient terms

    Assuming the liquids are incompressible, i.e., the speed of sound is infinite, the pressure-gradient term must be treated implicitly.
    See :ref:`the temporal integration of the momentum balance <momentum_integration>`.

********************************************
Approximate Factorization of Diffusive Terms
********************************************

Here we focus on the implicit treatment of the diffusive terms.
Applying :ref:`the combined scheme <time_marchers>` to :ref:`the momentum and internal energy balances <strong_conservation_form>` yields, for each Runge-Kutta sub-step:

.. math::

    \vel{i}^{k+1}
    -
    \vel{i}^{k}
    =
    \gamma^k
    \Delta t
    \frac{\sqrt{Pr}}{\sqrt{Ra}}
    \frac{1}{J}
    \dif{}{\gcs{j}}
    \left[
        \frac{J}{\sfact{j}}
        \frac{1}{\sfact{j}}
        \dif{}{\gcs{j}}
        \left\{
            c
            \vel{i}^{k+1}
            +
            \left( 1 - c \right)
            \vel{i}^k
        \right\}
    \right]
    +
    \left(
        \text{others}
    \right)_i,

and

.. math::

    T^{k+1}
    -
    T^{k}
    =
    \gamma^k
    \Delta t
    \frac{1}{\sqrt{Pr} \sqrt{Ra}}
    \frac{1}{J}
    \dif{}{\gcs{j}}
    \left[
        \frac{J}{\sfact{j}}
        \frac{1}{\sfact{j}}
        \dif{}{\gcs{j}}
        \left\{
            c
            T^{k+1}
            +
            \left( 1 - c \right)
            T^k
        \right\}
    \right]
    +
    \left(
        \text{others}
    \right),

respectively, where :math:`c` is a coefficient specifying the implicit treatment:

.. math::

    c
    =
    \begin{cases}
        \text{Euler explicit} & 0, \\
        \text{Crank-Nicolson} & \frac{1}{2}, \\
        \text{Euler implicit} & 1.
    \end{cases}

Here the last terms include the advective, pressure-gradient, and buoyancy terms which are not important here.
Since they are almost identical, we only focus on the temperature relation, which yields a Helmholtz equation:

.. math::

    \left\{
        1
        -
        c \gamma^k \Delta t
        \frac{1}{\sqrt{Pr} \sqrt{Ra}}
        \frac{1}{J}
        \dif{}{\gcs{j}}
        \left(
            \frac{J}{\sfact{j}}
            \frac{1}{\sfact{j}}
            \dif{}{\gcs{j}}
        \right)
    \right\}
    T^{k+1}
    =
    \left\{
        1
        +
        \left(
            1
            -
            c
        \right)
        \gamma^k \Delta t
        \frac{1}{\sqrt{Pr} \sqrt{Ra}}
        \frac{1}{J}
        \dif{}{\gcs{j}}
        \left(
            \frac{J}{\sfact{j}}
            \frac{1}{\sfact{j}}
            \dif{}{\gcs{j}}
        \right)
    \right\}
    T^{k}
    +
    \left(
        \text{others}
    \right).

Although this equation can be solved in a similar way as :ref:`solving Poisson equations <poisson_equation>`, we simplify it by utilising the approximate factorisation (|DUKOWICZ1992|) as follows.

First, we rewrite the equations as

.. math::

    \newcommand{\lap}[2]{
        {#2} \gamma^k \Delta t
        \frac{1}{\sqrt{Pr} \sqrt{Ra}}
        \frac{1}{J}
        \dif{}{\gcs{#1}}
        \left(
            \frac{J}{\sfact{#1}}
            \frac{1}{\sfact{#1}}
            \dif{}{\gcs{#1}}
        \right)
    }
    \left\{
        1
        -
        \lap{j}{c}
    \right\}
    \Delta T
    =
    \lap{j}{}
    T^k
    +
    \left(
        \text{others}
    \right),

where

.. math::

    \Delta T
    \equiv
    T^{n+1}
    -
    T^{n  }.

We approximate the left-hand side

.. math::

    \left\{
        1
        -
        \lap{1}{c}
        -
        \lap{2}{c}
        -
        \lap{3}{c}
    \right\}
    \Delta T

as

.. math::

    \left\{
        1
        -
        \lap{1}{c}
    \right\}
    \left\{
        1
        -
        \lap{2}{c}
    \right\}
    \left\{
        1
        -
        \lap{3}{c}
    \right\}
    \Delta T.

The leading-order error induced by this approximation is

.. math::

    \lap{1}{c}
    \Delta T
    \times
    \lap{2}{c}
    \Delta T
    =
    \mathcal{O} \left( \Delta t^3 \right)

by assuming

.. math::

    \Delta T
    \sim
    \Delta t.

Since :math:`\mathcal{O} \left( \Delta t^3 \right)` is comparable to the dominant error of :ref:`the three-step explicit Runge-Kutta scheme <time_marchers>`, this treatment is justified.
(Note that this is not accepted if we adopt a more accurate scheme to integrate the equation in time.)
Note that we need to solve :ref:`a linear system <linear_system>` here.

.. note::

    * Overhead

        In this project, the implicit treatment in the :math:`x` (wall-normal) direction can be easily achieved, whose cost is up to a few percent.
        The implicit treatments in the other directions, however, require certain amount of MPI communication, whose overhead can be more than :math:`100` percent.

    * Monotonicity

        Although the Crank-Nicolson scheme can eliminate the stability restriction, monotonicity (i.e., temperature is bounded between the two boundary values) is not guaranteed.
        Unfortunately, to guarantee the monotonicity, :math:`\Delta t \propto \left( \sfact{i} \right)^2` should be satisfied (see e.g., |HORVATH2000|).
        Although this restriction disappears by adopting the Euler implicit scheme (at the expense of the temporal accuracy), this issue is neglected in this project for now.

