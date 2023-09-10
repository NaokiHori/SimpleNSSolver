
.. include:: /references.txt

.. _time_marchers:

#####################
Time-marching schemes
#####################

*************
Problem setup
*************

I am interested in integrating several partial differential equations in this project.
When only the temporal derivatives are focused, I can write these equations as

.. math::

   \frac{df}{dt}
   =
   g

in general.

.. note::

   I assume :math:`g` is not a function of :math:`t` explicitly, i.e. the following discussion is limited for `autonomous systems <https://en.wikipedia.org/wiki/Autonomous_system_(mathematics)>`_.
   For instance,

   .. math::

      \der{T}{t}
      =
      -
      u_j \der{T}{x_j}
      +
      \frac{1}{\sqrt{Pr} \sqrt{Ra}} \der{}{x_j} \der{T}{x_j}

   includes no explicit :math:`t` dependency in the right-hand-side terms.

This is equal to

.. math::

   \int_{t^{n  }}^{t^{n+1}} \frac{df}{dt} dt
   =
   \int_{t^{n  }}^{t^{n+1}} g dt,

or

.. math::

   f^{n+1}
   -
   f^{n  }
   =
   \int_{t^{n  }}^{t^{n+1}} g dt,

where I introduce

.. math::

   f^{n  }
   \equiv
   f \left( t^{n  } \right).

The main focus of this page is how to approximate the right-hand side.
In other words, the objective is to approximate the Taylor-series expansion around :math:`t^{n  }`:

.. math::

   f^{n+1}
   &
   =
   f^{n  }
   +
   \frac{df^n}{dt} \Delta t
   +
   \frac{1}{2} \frac{d^2f^n}{dt^2} \left( \Delta t \right)^2
   +
   \mathcal{O} \left( \Delta t^3 \right) \\
   &
   =
   f^{n  }
   +
   g^{n  } \Delta t
   +
   \frac{1}{2} \left( \frac{dg}{dt} \right)^n \left( \Delta t \right)^2
   +
   \mathcal{O} \left( \Delta t^3 \right).

*******************
Zeroth-order scheme
*******************

The simplest way to approximate the above series expansion is

.. math::

   f^{n+1}
   \approx
   f^{n  },

where only the first term in the right-hand side is extracted.
The leading-order error (local truncation error) is

.. math::

   g^{n  } \Delta t
   =
   \mathcal{O} \left( \Delta t \right),

which indicates the first-order accuracy.

Since I am interested in studying the system after it has been integrated for a long time, I should focus on the global truncation error, which is :math:`\mathcal{O} \left( 1 \right)`.
This means that, even when I make :math:`\Delta t` smaller, the solution never converges and thus this scheme is useless.

.. note::

   Hereafter temporal accuracy is not based on the lobal truncation errors but the global ones.

********************
Euler-forward scheme
********************

A simple and practical scheme to integrate the above equation in time would be to extract up to the second term of the series expansion:

.. math::

   f^{n+1}
   =
   f^{n  }
   +
   g^{n  }
   \Delta t,

or equivalently

.. math::

   \Delta f
   & =
   g^{n  }
   \Delta t, \\
   f^{n+1}
   & =
   f^{n  }
   +
   \Delta f,

which is the Euler-forward (or Euler-explicit) scheme.

Since the deviation from the series expansion is :math:`\mathcal{O} \left( \Delta t^2 \right)`, this scheme has the first-order accuracy in time.

Although it is very simple and easy to understand, there are mainly two issues in the above scheme:

   #. It only has the first-order accuracy in time.

   #. All terms are treated explicitly.

They are elaborated in the following sections.

*********************************
Fully-explicit Runge-Kutta scheme
*********************************

To make the temporal accuracy second-order, I need to take three terms in the series expansion above.
In this project, I adopt the explicit Runge-Kutta scheme.
Since the system is autonomous, I use a three-step and low-storage scheme to achieve the second-order accuracy in time:

.. math::

   &\text{do}\,\,\,\, k = 0, 2 \\
   &\,\,\,\,
   \Delta f
   =
   \left(
      \alpha^k g^{k  }
      +
      \beta^k  g^{k-1}
   \right) \Delta t, \\
   &\,\,\,\,
   f^{k+1}
   =
   f^{k  }
   +
   \Delta f. \\
   &\text{enddo}

Although :math:`\alpha^k` and :math:`\beta^k`, which are coefficients to achieve the second-order accuracy in time, `have several possible combinations <https://en.wikipedia.org/wiki/List_of_Runge–Kutta_methods>`_, in this project, I adopt

.. math::

   \left(\alpha^0, \alpha^1, \alpha^2 \right)
   =
   \left(32/60, 25/60, 45/60 \right),

and

.. math::

   \left(\beta^0, \beta^1, \beta^2 \right)
   =
   \left(0, -17/60, -25/60 \right),

following e.g. |RAI1991|, |VERZICCO1996|, |COSTA2018|.

Although not used in the above equation, I introduce

.. math::

   \gamma^{k  }
   \equiv
   \alpha^{k  }
   +
   \beta^{k  },

or explicitly

.. math::

   \left(\gamma^0, \gamma^1, \gamma^2 \right)
   =
   \left(32/60, 8/60, 20/60 \right),

which is used in the next section.

.. note::

   An important takeaway here is that, since

   .. math::

      \sum_{k = 0}^{2} \gamma^k
      =
      1,

   I can regard the above three-step Runge-Kutta scheme as the combination of the three Euler-forward scheme having different and smaller time steps :math:`\gamma^k \Delta t`.

.. mydetails:: Proof of the second-order accuracy

   Because

   .. math::

      f^{k+1}
      =
      f^{k  }
      +
      \alpha^k g^{k  } \Delta t
      +
      \beta^k  g^{k-1} \Delta t,

   I find

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

   By using this relation repeatedly, I have

   .. math::

      f^1 &= f^n + \alpha^0 g^n \Delta t \\
          &= f^n + \alpha^0 \frac{df^n}{dt} \Delta t, \\
      g^1 &= \frac{df^n}{dt} + \alpha^0 \frac{d^2f^n}{dt^2} \Delta t, \\
      f^2 &= f^1 + \alpha^1 g^1 \Delta t + \beta^1 g^0 \Delta t \\
          &= f^n + \alpha^0 \frac{df^n}{dt} \Delta t + \alpha^1 \left( \frac{df^n}{dt} + \alpha^0 \frac{d^2f^n}{dt^2} \Delta t \right) \Delta t + \beta^1 \frac{df^n}{dt} \Delta t, \\
          &= f^n + \left( \alpha^0 + \alpha^1 + \beta^1 \right) \frac{df^n}{dt} \Delta t + \alpha^0 \alpha^1 \frac{d^2f^n}{dt^2} \Delta t^2, \\
      g^2 &= \frac{df^n}{dt} + \left( \alpha^0 + \alpha^1 + \beta^1 \right) \frac{d^2f^n}{dt^2} \Delta t + \alpha^0 \alpha^1 \frac{d^3f^n}{dt^3} \Delta t^2, \\
      f^3 &= f^{n+1} \\
          &= f^2 + \alpha^2 g^2 \Delta t + \beta^2 g^1 \Delta t \\
          &= f^n + \left( \alpha^0 + \alpha^1 + \beta^1 \right) \frac{df^n}{dt} \Delta t + \alpha^0 \alpha^1 \frac{d^2f^n}{dt^2} \Delta t^2 \\
          &+ \alpha^2 \left( \frac{df^n}{dt} + \left( \alpha^0 + \alpha^1 + \beta^1 \right) \frac{d^2f^n}{dt^2} \Delta t + \alpha^0 \alpha^1 \frac{d^3f^n}{dt^3} \Delta t^2 \right) \Delta t \\
          &+ \beta^2 \left( \frac{df^n}{dt} + \alpha^0 \frac{d^2f^n}{dt^2} \Delta t \right) \Delta t, \\
          &= f^n + \left( \alpha^0 + \alpha^1 + \beta^1 + \alpha^2 + \beta^2 \right) \frac{df^n}{dt} \Delta t \\
          &+ \left( \alpha^0 \alpha^1 + \alpha^0 \alpha^2 + \alpha^1 \alpha^2 + \alpha^2 \beta^1 + \alpha^0 \beta^2 \right) \frac{d^2f^n}{dt^2} \Delta t^2 \\
          &+ \mathcal{O} \left( \Delta t^3 \right) \\
          &= f^n + \frac{df^n}{dt} \Delta t + \frac{1}{2} \frac{d^2f^n}{dt^2} \Delta t^2 + \mathcal{O} \left( \Delta t^3 \right),

   i.e.

   .. math::

      f^{n+1}
      =
      f^n
      +
      \frac{df^n}{dt} \Delta t
      +
      \frac{1}{2} \frac{d^2f^n}{dt^2} \left( \Delta t^2 \right)
      +
      \mathcal{O} \left( \Delta t^3 \right),

   and thus this scheme has the second-order accuracy in time.

***********************
Semi-implicit treatment
***********************

As discussed in :ref:`the implicit treatment of the diffusive terms <implicit_treatment>`, I often need to treat some terms implicitly in time, i.e.

.. math::

   \Delta f
   \equiv
   f^{n+1}
   -
   f^{n  }
   =
   \frac{1}{2}
   g^{n  }
   \Delta t
   +
   \frac{1}{2}
   g^{n+1}
   \Delta t,

where `Crank-Nicolson scheme <https://en.wikipedia.org/wiki/Crank–Nicolson_method>`_ is adopted to keep the second-order accuracy in time.

.. mydetails:: Proof of the second-order accuracy

   By using

   .. math::

      g^{n+1}
      = \frac{df^{n+1}}{dt}
      = \frac{df^n}{dt} + \frac{d^2f^n}{dt^2} \Delta t + \frac{1}{2} \frac{d^3f^n}{dt^3} \Delta t^2 + \mathcal{O} \left( \Delta t^3 \right),

   I find

   .. math::

      f^{n+1}
      &= f^n + \frac{1}{2} \left( \frac{df^n}{dt} + \frac{d^2f^n}{dt^2} \Delta t + \frac{1}{2} \frac{d^3f^n}{dt^3} \Delta t^2 + \frac{df^n}{dt} \right) \Delta t \\
      &= f^n + \frac{df^n}{dt} \Delta t + \frac{1}{2} \frac{d^2f^n}{dt^2} \Delta t^2 + \mathcal{O} \left( \Delta t^3 \right),

   i.e. second-order accuracy in time.

To embed this in the above Runge-Kutta iterations, I regard the above three-step Runge-Kutta method as the combination of three Euler-forward iterations with smaller time step sizes, i.e.

.. math::

   f^{1}
   -
   f^{0}
   =
   \left( \gamma^0 \Delta t \right) g^{0},

.. math::

   f^{2}
   -
   f^{1}
   =
   \left( \gamma^1 \Delta t \right) g^{1},

.. math::

   f^{3}
   -
   f^{2}
   =
   \left( \gamma^2 \Delta t \right) g^{2},

since

.. math::

   \sum_{k = 0}^{2} \gamma^k
   =
   1.

Then I replace the Euler-forward scheme with the Crank-Nicolson scheme, giving

.. math::

   &\text{do}\,\,\,\, k = 0, 2 \\
   &\,\,\,\,
   \Delta f
   =
   \frac{1}{2}
   \left(
      g^{k  }
      +
      g^{k+1}
   \right) \gamma^k \Delta t, \\
   &\,\,\,\,
   f^{k+1}
   =
   f^{k  }
   +
   \Delta f. \\
   &\text{enddo}

In summary, the conclusive scheme leads to

.. math::

   &\text{do}\,\,\,\, k = 0, 2 \\
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

