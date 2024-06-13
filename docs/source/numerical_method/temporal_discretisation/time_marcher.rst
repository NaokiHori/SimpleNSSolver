
.. _time_marchers:

.. include:: /references.txt

#####################
Time-Marching Schemes
#####################

*************
Problem Setup
*************

We focus on integrating several partial differential equations involving both temporal and spatial derivatives in this project.
When focusing only on the temporal derivatives, all equations can be written as

.. math::

    \newcommand{\oldvar}[1]{{#1^{n}}}
    \newcommand{\newvar}[1]{{#1^{n+1}}}
    \newcommand{\oerror}[1]{\mathcal{O} \left( \Delta t^{#1} \right)}
    \frac{df}{dt}
    =
    g.

.. note::

    We assume the right-hand-side term :math:`g` is not a function of time explicitly, i.e., the following discussion is limited for `autonomous systems <https://en.wikipedia.org/wiki/Autonomous_system_(mathematics)>`_.
    This applies to all equations discussed in this project.

This is equivalent to

.. math::

    \int_\oldvar{t}^\newvar{t} \frac{df}{dt} dt
    =
    \int_\oldvar{t}^\newvar{t} g dt,

or

.. math::

    \newvar{f}
    -
    \oldvar{f}
    =
    \int_\oldvar{t}^\newvar{t} g dt,

where we introduce

.. math::

    \oldvar{f}
    &
    \equiv
    f \left( t = \oldvar{t} \right),

    \newvar{f}
    &
    \equiv
    f \left( t = \newvar{t} \right),

and

.. math::

    \Delta t
    \equiv
    \newvar{t}
    -
    \oldvar{t}.

To obtain :math:`\newvar{f}` from :math:`\oldvar{f}`, we need to evaluate the right-hand-side integral.
In this page, we aim at approximating the right-hand-side term using the Taylor-series expansion around :math:`\oldvar{t}`:

.. math::

    \newvar{f}
    &
    =
    \oldvar{f}
    +
    \frac{d\oldvar{f}}{dt} \Delta t
    +
    \frac{1}{2} \frac{d^2\oldvar{f}}{dt^2} \left( \Delta t \right)^2
    +
    \oerror{3}

    &
    =
    \oldvar{f}
    +
    \oldvar{g} \Delta t
    +
    \frac{1}{2} \frac{d\oldvar{g}}{dt} \left( \Delta t \right)^2
    +
    \oerror{3}.

*******************
Zeroth-Order Scheme
*******************

The simplest but useless way to approximate the above Taylor series expansion is

.. math::

    \newvar{f}
    \approx
    \oldvar{f},

where we only pick-up the first term.
The leading-order error (local truncation error) is

.. math::

    \oldvar{g} \Delta t
    =
    \oerror{1},

which indicates the first-order accuracy.
Since our objective is integrate the equation for a long time, we should consider the global truncation error, giving a reduced order of accuracy: :math:`\oerror{0}`.
This indicates that the error never shrinks even with a smaller :math:`\Delta t`, which is useless.

.. note::

    Hereafter we only focus on the global truncation errors.

********************
Euler-Forward Scheme
********************

We obtain a simple but (to some extent) practical scheme by extracting the first two terms:

.. math::

    \newvar{f}
    =
    \oldvar{f}
    +
    \oldvar{g}
    \Delta t,

or equivalently

.. math::

    \Delta f
    &
    =
    \oldvar{g}
    \Delta t,

    \newvar{f}
    &
    =
    \oldvar{f}
    +
    \Delta f,

which is known as the Euler-forward (or Euler-explicit) scheme.

Since the deviation from the series expansion is :math:`\oerror{2}` (locally), this scheme has the first-order accuracy (globally).

Although it is very simple and easy to use, two clear issues exist:

   * It only has the first-order accuracy in time,

   * All terms are treated explicitly.

*********************************
Fully-Explicit Runge-Kutta Scheme
*********************************

To achieve a higher-order accuracy, we need to extract more than two terms from the series expansion above.
In this project, we adopt the explicit Runge-Kutta scheme.
Since the system is autonomous, we use a three-step and low-storage scheme to achieve the second-order accuracy in time:

.. math::

    &
    \text{do}\,\,k = 0, 2

    &
    \,\,\,\,
    \Delta f
    =
    \left(
       \alpha^k g^{k  }
       +
       \beta^k  g^{k-1}
    \right) \Delta t,

    &
    \,\,\,\,
    f^{k+1}
    =
    f^{k  }
    +
    \Delta f.

    &
    \text{enddo}

Although the coefficients :math:`\alpha^k,\beta^k` have `several possibilities <https://en.wikipedia.org/wiki/List_of_Runge–Kutta_methods>`_, we adopt

.. math::

   \left( \alpha^0, \alpha^1, \alpha^2 \right)
   =
   \left( \frac{32}{60}, \frac{25}{60}, \frac{45}{60} \right),

and

.. math::

   \left( \beta^0, \beta^1, \beta^2 \right)
   =
   \left( \frac{0}{60}, -\frac{17}{60}, -\frac{25}{60} \right),

following e.g., |RAI1991|, |VERZICCO1996|, |COSTA2018|.
For later convenience, we also introduce

.. math::

   \gamma^{k  }
   \equiv
   \alpha^{k  }
   +
   \beta^{k  },

or explicitly

.. math::

   \left( \gamma^0, \gamma^1, \gamma^2 \right)
   =
   \left( \frac{32}{60}, \frac{8}{60}, \frac{20}{60} \right).

Note that, since we have

.. math::

  \sum_{k = 0}^{2} \gamma^k
  =
  1,

the above three-step Runge-Kutta scheme is essentially a combination of three Euler-forward scheme with :math:`\gamma^k \Delta t` as time-step sizes.

.. mydetails:: Proof of the second-order accuracy

    .. include:: derivation/rk.rst

******************
Implicit Treatment
******************

As justified :ref:`later <implicit_treatment>`, we treat some terms implicitly in time; among others we adopt `Crank-Nicolson scheme <https://en.wikipedia.org/wiki/Crank–Nicolson_method>`_:

.. math::

    \Delta f
    &
    =
    \frac{1}{2}
    \oldvar{g}
    \Delta t
    +
    \frac{1}{2}
    g^{n+1}
    \Delta t,

    \newvar{f}
    -
    \oldvar{f}
    &
    =
    \Delta f,

which has the second-order accuracy in time as well.

.. mydetails:: Proof of the second-order accuracy

    .. include:: derivation/cn.rst

***************
Combined scheme
***************

We embed the implicit treatment (Crank-Nicolson scheme) into the Runge-Kutta iterations:

.. math::

    &\text{do}\,\, k = 0, 2 \\
    &\,\,\,\,
    \Delta f
    =
    \left(
      \alpha^k g^{k  }
      +
      \beta^k  g^{k-1}
    \right) \Delta t
    +
    \frac{1}{2}
    \left(
      h^{k  }
      +
      h^{k+1}
    \right) \gamma^k \Delta t, \\
    &\,\,\,\,
    f^{k+1}
    =
    f^{k  }
    +
    \Delta f, \\
    &\text{enddo}

where :math:`g` and :math:`h` are terms treated explicitly and implicitly in time, respectively.

