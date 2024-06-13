Because of

.. math::

    g^{n+1}
    =
    \frac{d\newvar{f}}{dt}
    =
    \frac{df^n}{dt} + \frac{d^2f^n}{dt^2} \Delta t + \frac{1}{2} \frac{d^3f^n}{dt^3} \Delta t^2 + \oerror{3},

we find

.. math::

    \newvar{f}
    &
    =
    f^n
    +
    \frac{1}{2}
    \frac{df^n}{dt}
    \Delta t
    +
    \frac{1}{2} \left( \frac{df^n}{dt} + \frac{d^2f^n}{dt^2} \Delta t + \frac{1}{2} \frac{d^3f^n}{dt^3} \Delta t^2 + \frac{df^n}{dt} \right) \Delta t

    &
    =
    f^n + \frac{df^n}{dt} \Delta t + \frac{1}{2} \frac{d^2f^n}{dt^2} \Delta t^2 + \oerror{3},

indicating that this scheme has the third-order accuracy (locally) and the second-order accuracy (globally).

