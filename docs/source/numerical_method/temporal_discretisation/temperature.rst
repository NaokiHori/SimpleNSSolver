
.. _temperature_integration:

.. include:: /references.txt

#######################
Internal Energy Balance
#######################

The equation of the internal energy is

.. math::

    \pder{T}{t}
    =
    A
    +
    D_x
    +
    D_y
    +
    D_z,

where :math:`A` is the advective terms:

.. math::

    A
    \equiv
    -
    \dtempadv{1}
    -
    \dtempadv{2}
    -
    \dtempadv{3},

while :math:`D_i` is the diffusive term involving spatial differentiation in the :math:`i`-th direction:

.. math::

    D_i
    \equiv
    \dtempdif{i}
    \,\,
    (\text{no summation over}\,i).

See :ref:`the spatial discretization <discrete_internal_energy_balance>`.

The temporal discretization for each Runge-Kutta iteration leads to

.. math::

    \Delta T
    &
    =
    \alpha^k \Delta t \left( A^{k  } + D^{k  } \right)
    +
    \beta^k  \Delta t \left( A^{k-1} + D^{k-1} \right),

    T^{k+1}
    &
    =
    T^k
    +
    \Delta T,

when all diffusive terms are treated explicitly, while

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
    \Delta T
    &
    =
    \alpha^k \Delta t A^{k  }
    +
    \beta^k  \Delta t A^{k-1}
    +
    \gamma^k \Delta t D^{k  },

    T^{k+1}
    &
    =
    T^k
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
    \Delta T,

when all diffusive terms are treated implicitly.

Diffusive terms are sometimes partially treated implicitly (c.f., |VANDERPOEL2015|).
When only the diffusive term in :math:`x` direction is implicitly treated (the default configuration in this project), we have

.. math::

    \Delta T
    &
    =
    \alpha^k \Delta t \left( A^{k  } + D_2^{k  } + D_3^{k  } \right)
    +
    \beta^k  \Delta t \left( A^{k-1} + D_2^{k-1} + D_3^{k-1} \right)
    +
    \gamma^k \Delta t D_1^{k  },

    T^{k+1}
    &
    =
    T^k
    +
    \left\{
        1
        -
        \lap{1}{c}
    \right\}^{-1}
    \Delta T.

First, the explicit and implicit terms are calculated and stored to the corresponding buffers:

.. myliteralinclude:: /../../src/fluid/predict/t.c
    :language: c
    :tag: compute right-hand-side terms, which are added to buffers

The buffers are used to compute :math:`\Delta T`:

.. myliteralinclude:: /../../src/fluid/predict/t.c
    :language: c
    :tag: compute increments

When necessary, linear systems are solved to take care of the implicit treatments:

.. myliteralinclude:: /../../src/fluid/predict/t.c
    :language: c
    :tag: solve linear systems in x

.. myliteralinclude:: /../../src/fluid/predict/t.c
    :language: c
    :tag: solve linear systems in y

.. myliteralinclude:: /../../src/fluid/predict/t.c
    :language: c
    :tag: solve linear systems in z

Finally the temperature field is updated:

.. myliteralinclude:: /../../src/fluid/predict/t.c
    :language: c
    :tag: update temperature field

