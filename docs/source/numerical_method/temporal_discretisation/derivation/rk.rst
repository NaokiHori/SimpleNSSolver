Because of

.. math::

    f^{k+1}
    =
    f^{k  }
    +
    \alpha^k g^{k  } \Delta t
    +
    \beta^k  g^{k-1} \Delta t,

we find

.. math::

    g^{k+1}
    =
    \frac{df^{k+1}}{dt}
    =
    \frac{df^k}{dt}
    +
    \alpha^k \frac{d^2f^k}{dt^2} \Delta t
    +
    \beta^k \frac{d^2f^{k-1}}{dt^2} \Delta t.

We use this relation repeatedly.
For :math:`k = 1`, we have

.. math::

    f^1
    =
    f^n
    +
    \alpha^0 g^n \Delta t
    =
    f^n + \alpha^0 \frac{df^n}{dt} \Delta t,

whose derivation leads to

.. math::

    g^1
    =
    \frac{df^n}{dt}
    +
    \alpha^0 \frac{d^2f^n}{dt^2} \Delta t.

Note that :math:`k = 0` corresponds to :math:`n` (old information).
For :math:`k = 2`, we have

.. math::

    f^2
    &
    =
    f^1
    +
    \alpha^1 g^1 \Delta t + \beta^1 g^0 \Delta t

    &
    =
    f^n
    +
    \alpha^0 \frac{df^n}{dt} \Delta t + \alpha^1 \left( \frac{df^n}{dt} + \alpha^0 \frac{d^2f^n}{dt^2} \Delta t \right) \Delta t + \beta^1 \frac{df^n}{dt} \Delta t,

    &
    =
    f^n + \left( \alpha^0 + \alpha^1 + \beta^1 \right) \frac{df^n}{dt} \Delta t + \alpha^0 \alpha^1 \frac{d^2f^n}{dt^2} \Delta t^2,

whose derivation leads to

.. math::

    g^2
    =
    \frac{df^n}{dt} + \left( \alpha^0 + \alpha^1 + \beta^1 \right) \frac{d^2f^n}{dt^2} \Delta t + \alpha^0 \alpha^1 \frac{d^3f^n}{dt^3} \Delta t^2.

For :math:`k = 3`, we have

.. math::

    f^3
    &
    =
    f^2 + \alpha^2 g^2 \Delta t + \beta^2 g^1 \Delta t

    &
    =
    f^n + \left( \alpha^0 + \alpha^1 + \beta^1 \right) \frac{df^n}{dt} \Delta t + \alpha^0 \alpha^1 \frac{d^2f^n}{dt^2} \Delta t^2
    +
    \alpha^2 \left( \frac{df^n}{dt} + \left( \alpha^0 + \alpha^1 + \beta^1 \right) \frac{d^2f^n}{dt^2} \Delta t + \alpha^0 \alpha^1 \frac{d^3f^n}{dt^3} \Delta t^2 \right) \Delta t
    +
    \beta^2 \left( \frac{df^n}{dt} + \alpha^0 \frac{d^2f^n}{dt^2} \Delta t \right) \Delta t,

    &
    =
    f^n + \left( \alpha^0 + \alpha^1 + \beta^1 + \alpha^2 + \beta^2 \right) \frac{df^n}{dt} \Delta t
    +
    \left( \alpha^0 \alpha^1 + \alpha^0 \alpha^2 + \alpha^1 \alpha^2 + \alpha^2 \beta^1 + \alpha^0 \beta^2 \right) \frac{d^2f^n}{dt^2} \Delta t^2
    +
    \oerror{3}

    &
    =
    f^n + \frac{df^n}{dt} \Delta t + \frac{1}{2} \frac{d^2f^n}{dt^2} \Delta t^2 + \oerror{3}.

Since :math:`k = 3` corresponds to :math:`n + 1` (new information), we find

.. math::

    \newvar{f}
    =
    f^n
    +
    \frac{df^n}{dt} \Delta t
    +
    \frac{1}{2} \frac{d^2f^n}{dt^2} \left( \Delta t^2 \right)
    +
    \oerror{3},

indicating that this scheme has the third-order accuracy (locally) and the second-order accuracy (globally).

