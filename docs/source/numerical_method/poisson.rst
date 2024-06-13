
.. _poisson_equation:

################
Poisson Equation
################

To :ref:`integrate a momentum field in time <momentum_integration>` (specifically to enforce :ref:`the incompressibility <discrete_incompressibility>`), we need to solve a Poisson equation:

.. math::

    \frac{1}{J}
    \dif{}{\gcs{1}}
    \left(
        \frac{J}{\sfact{1}}
        \frac{1}{\sfact{1}}
        \dif{p}{\gcs{1}}
    \right)
    +
    \frac{1}{J}
    \dif{}{\gcs{2}}
    \left(
        \frac{J}{\sfact{2}}
        \frac{1}{\sfact{2}}
        \dif{p}{\gcs{2}}
    \right)
    +
    \frac{1}{J}
    \dif{}{\gcs{3}}
    \left(
        \frac{J}{\sfact{3}}
        \frac{1}{\sfact{3}}
        \dif{p}{\gcs{3}}
    \right)
    =
    q,

where :math:`p` and :math:`q` are both defined at cell centers :math:`\left( \ccidx{i}, \ccidx{j}, \ccidx{k} \right)`.

Since directions except :math:`x` are homogeneous (i.e., periodic boundary conditions are imposed and the grid sizes are equal), we consider the spectral representation (discrete backward Fourier transform) of :math:`p`:

.. math::

    \newcommand{ \wavy}{\exp \left(   I \frac{2 \pi}{\ngp{2}} j m        \right)}
    \newcommand{ \wavz}{\exp \left(   I \frac{2 \pi}{\ngp{3}} k n        \right)}
    \newcommand{\iwavy}{\exp \left( - I \frac{2 \pi}{\ngp{2}} j m^\prime \right)}
    \newcommand{\iwavz}{\exp \left( - I \frac{2 \pi}{\ngp{3}} k n^\prime \right)}
    \vat{p}{\ccidx{i}, \ccidx{j}, \ccidx{k}}
    =
    \sum_{n = 0}^{\ngp{3} - 1}
    \sum_{m = 0}^{\ngp{2} - 1}
    \vat{P}{\ccidx{i}, \ccidx{m}, \ccidx{n}}
    \wavy
    \wavz,

where :math:`I` is the imaginary unit :math:`\sqrt{-1}`.
Note that the discrete forward Fourier transform yields

.. math::

    \sum_{k = 0}^{\ngp{3} - 1}
    \sum_{j = 0}^{\ngp{2} - 1}
    \vat{p}{\ccidx{i}, \ccidx{j}, \ccidx{k}}
    \iwavy
    \iwavz
    &
    =
    \sum_{k = 0}^{\ngp{3} - 1}
    \sum_{j = 0}^{\ngp{2} - 1}
    \sum_{n = 0}^{\ngp{3} - 1}
    \sum_{m = 0}^{\ngp{2} - 1}
    \vat{P}{\ccidx{i}, \ccidx{m}, \ccidx{n}}
    \wavy
    \wavz
    \iwavy
    \iwavz

    &
    =
    \sum_{k = 0}^{\ngp{3} - 1}
    \sum_{j = 0}^{\ngp{2} - 1}
    \sum_{n = 0}^{\ngp{3} - 1}
    \sum_{m = 0}^{\ngp{2} - 1}
    \vat{P}{\ccidx{i}, \ccidx{m}, \ccidx{n}}
    \delta_{m m^\prime}
    \delta_{n n^\prime}

    &
    =
    \ngp{2}
    \ngp{3}
    \sum_{n = 0}^{\ngp{3} - 1}
    \sum_{m = 0}^{\ngp{2} - 1}
    \vat{P}{\ccidx{i}, \ccidx{m}, \ccidx{n}}
    \delta_{m m^\prime}
    \delta_{n n^\prime}

    &
    =
    \ngp{2}
    \ngp{3}
    \vat{P}{\ccidx{i}, \ccidx{m^\prime}, \ccidx{n^\prime}}.

Since we have

.. math::

    \vat{p}{\ccidx{i} \pm 1, \ccidx{j}, \ccidx{k}}
    &
    =
    \sum_{n = 0}^{\ngp{3} - 1}
    \sum_{m = 0}^{\ngp{2} - 1}
    \vat{P}{\ccidx{i} \pm 1, \ccidx{m}, \ccidx{n}}
    \wavy
    \wavz,

    \vat{p}{\ccidx{i}, \ccidx{j} \pm 1, \ccidx{k}}
    &
    =
    \sum_{n = 0}^{\ngp{3} - 1}
    \sum_{m = 0}^{\ngp{2} - 1}
    \vat{P}{\ccidx{i}, \ccidx{m}, \ccidx{n}}
    \exp \left( \pm I \frac{2 \pi}{\ngp{2}} m \right)
    \wavy
    \wavz,

    \vat{p}{\ccidx{i}, \ccidx{j}, \ccidx{k} \pm 1}
    &
    =
    \sum_{n = 0}^{\ngp{3} - 1}
    \sum_{m = 0}^{\ngp{2} - 1}
    \vat{P}{\ccidx{i}, \ccidx{m}, \ccidx{n}}
    \exp \left( \pm I \frac{2 \pi}{\ngp{3}} n \right)
    \wavy
    \wavz,

the second-order derivatives in the Poisson equation can be reformulated as

.. math::

    \newcommand{\coefl}{C_l}
    \newcommand{\coefu}{C_u}
    \frac{1}{J}
    \dif{}{\gcs{1}}
    \left(
        \frac{J}{\sfact{1}}
        \frac{1}{\sfact{1}}
        \dif{p}{\gcs{1}}
    \right)
    &
    =
    \coefl
    \vat{p}{\ccidx{i} - 1, \ccidx{j}, \ccidx{k}}
    -
    \left(
        \coefl
        +
        \coefu
    \right)
    \vat{p}{\ccidx{i}, \ccidx{j}, \ccidx{k}}
    +
    \coefu
    \vat{p}{\ccidx{i} + 1, \ccidx{j}, \ccidx{k}}

    &
    =
    \coefl
    \sum_{n = 0}^{\ngp{3} - 1}
    \sum_{m = 0}^{\ngp{2} - 1}
    \vat{P}{\ccidx{i} - 1, \ccidx{m}, \ccidx{n}}
    \wavy
    \wavz
    -
    \left(
        \coefl
        +
        \coefu
    \right)
    \sum_{n = 0}^{\ngp{3} - 1}
    \sum_{m = 0}^{\ngp{2} - 1}
    \vat{P}{\ccidx{i}, \ccidx{m}, \ccidx{n}}
    \wavy
    \wavz
    +
    \coefu
    \sum_{n = 0}^{\ngp{3} - 1}
    \sum_{m = 0}^{\ngp{2} - 1}
    \vat{P}{\ccidx{i} + 1, \ccidx{m}, \ccidx{n}}
    \wavy
    \wavz,

    \frac{1}{J}
    \dif{}{\gcs{2}}
    \left(
        \frac{J}{\sfact{2}}
        \frac{1}{\sfact{2}}
        \dif{p}{\gcs{2}}
    \right)
    &
    =
    \frac{1}{\sfact{2}}
    \frac{1}{\sfact{2}}
    \left(
        \vat{p}{\ccidx{i}, \ccidx{j} - 1, \ccidx{k}}
        -
        2
        \vat{p}{\ccidx{i}, \ccidx{j}    , \ccidx{k}}
        +
        \vat{p}{\ccidx{i}, \ccidx{j} + 1, \ccidx{k}}
    \right)

    &
    =
    \frac{1}{\sfact{2}}
    \frac{1}{\sfact{2}}
    \sum_{n = 0}^{\ngp{3} - 1}
    \sum_{m = 0}^{\ngp{2} - 1}
    \vat{P}{\ccidx{i}, \ccidx{m}, \ccidx{n}}
    \wavy
    \wavz
    \left\{
        \exp \left(   I \frac{2 \pi}{\ngp{2}} m \right)
        -
        2
        +
        \exp \left( - I \frac{2 \pi}{\ngp{2}} m \right)
    \right\}

    &
    =
    -
    \frac{1}{\sfact{2}}
    \frac{1}{\sfact{2}}
    \sum_{n = 0}^{\ngp{3} - 1}
    \sum_{m = 0}^{\ngp{2} - 1}
    4
    \vat{P}{\ccidx{i}, \ccidx{m}, \ccidx{n}}
    \wavy
    \wavz
    \sin^2 \left( \frac{\pi}{\ngp{2}} m \right),

    \frac{1}{J}
    \dif{}{\gcs{3}}
    \left(
        \frac{J}{\sfact{3}}
        \frac{1}{\sfact{3}}
        \dif{p}{\gcs{3}}
    \right)
    &
    =
    -
    \frac{1}{\sfact{3}}
    \frac{1}{\sfact{3}}
    \sum_{n = 0}^{\ngp{3} - 1}
    \sum_{m = 0}^{\ngp{2} - 1}
    4
    \vat{P}{\ccidx{i}, \ccidx{m}, \ccidx{n}}
    \wavy
    \wavz
    \sin^2 \left( \frac{\pi}{\ngp{3}} n \right),

where

.. math::

    \coefl
    &
    \equiv
    \frac{
        1
    }{
        \vat{J}{\ccidx{i}}
    }
    \frac{
        \vat{J}{\cmidx{i}}
    }{
        \vat{\sfact{1}}{\cmidx{i}}
    }
    \frac{
        1
    }{
        \vat{\sfact{1}}{\cmidx{i}}
    },

    \coefu
    &
    \equiv
    \frac{
        1
    }{
        \vat{J}{\ccidx{i}}
    }
    \frac{
        \vat{J}{\cpidx{i}}
    }{
        \vat{\sfact{1}}{\cpidx{i}}
    }
    \frac{
        1
    }{
        \vat{\sfact{1}}{\cpidx{i}}
    },

which are computed here:

.. myliteralinclude:: /../../src/fluid/compute_potential.c
    :language: c
    :tag: initialise tri-diagonal matrix in x direction

Consequently, we obtain

.. math::

    &
    \sum_{n = 0}^{\ngp{3} - 1}
    \sum_{m = 0}^{\ngp{2} - 1}
    \coefl
    \vat{P}{\ccidx{i} - 1, \ccidx{m}, \ccidx{n}}
    \wavy
    \wavz

    +
    &
    \sum_{n = 0}^{\ngp{3} - 1}
    \sum_{m = 0}^{\ngp{2} - 1}
    \left\{
        -
        \left(
            \coefl
            +
            \coefu
        \right)
        -
        4
        \frac{1}{\sfact{2}}
        \frac{1}{\sfact{2}}
        \sin^2 \left( \frac{\pi}{\ngp{2}} m \right)
        -
        4
        \frac{1}{\sfact{3}}
        \frac{1}{\sfact{3}}
        \sin^2 \left( \frac{\pi}{\ngp{3}} n \right)
    \right\}
    \vat{P}{\ccidx{i}, \ccidx{m}, \ccidx{n}}
    \wavy
    \wavz

    +
    &
    \sum_{n = 0}^{\ngp{3} - 1}
    \sum_{m = 0}^{\ngp{2} - 1}
    \coefu
    \vat{P}{\ccidx{i} + 1, \ccidx{m}, \ccidx{n}}
    \wavy
    \wavz

    =
    &
    \vat{q}{\ccidx{i}, \ccidx{j}, \ccidx{k}}.

By applying the discrete forward Fourier transform, we obtain

.. math::

    &
    \coefl
    \vat{P}{\ccidx{i} - 1, \ccidx{m}, \ccidx{n}}

    +
    &
    \left\{
        -
        \left(
            \coefl
            +
            \coefu
        \right)
        -
        4
        \frac{1}{\sfact{2}}
        \frac{1}{\sfact{2}}
        \sin^2 \left( \frac{\pi}{\ngp{2}} m \right)
        -
        4
        \frac{1}{\sfact{3}}
        \frac{1}{\sfact{3}}
        \sin^2 \left( \frac{\pi}{\ngp{3}} n \right)
    \right\}
    \vat{P}{\ccidx{i}, \ccidx{m}, \ccidx{n}}

    +
    &
    \coefu
    \vat{P}{\ccidx{i} + 1, \ccidx{m}, \ccidx{n}}

    =
    &
    \frac{1}{\ngp{2}}
    \frac{1}{\ngp{3}}
    \sum_{k = 0}^{\ngp{3} - 1}
    \sum_{j = 0}^{\ngp{2} - 1}
    \vat{q}{\ccidx{i}, \ccidx{j}, \ccidx{k}}
    \exp \left( - I \frac{2 \pi}{\ngp{2}} j m \right)
    \exp \left( - I \frac{2 \pi}{\ngp{3}} k n \right).

The resulting :ref:`tri-diagonal matrix <tdm>` for each :math:`m` and :math:`n` is solved here:

.. myliteralinclude:: /../../src/fluid/compute_potential.c
    :language: c
    :tag: solve tri-diagonal matrices

The sinusoidal functions (wave numbers) are computed here:

.. myliteralinclude:: /../../src/fluid/compute_potential.c
    :language: c
    :tag: y eigenvalues

.. myliteralinclude:: /../../src/fluid/compute_potential.c
    :language: c
    :tag: z eigenvalues

The overall procedure can be found here:

.. myliteralinclude:: /../../src/fluid/compute_potential.c
    :language: c
    :tag: solve Poisson equation

