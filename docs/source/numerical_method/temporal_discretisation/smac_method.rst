
.. include:: /references.txt

.. _smac_method:

##################################################
Integrating mass conservation and momentum balance
##################################################

The momentum equation in non-dimensional form is

.. math::

   \der{u_i}{t}
   =
   P_i
   +
   A_i
   +
   D_i,

where I introduce some symbols for notational convenience:

.. math::

   P_i^k
   \equiv
   -
   \dder{p^k}{x_i},

.. math::

   A_i^k
   \equiv
   -
   u_j^k \dder{u_i^k}{x_j},

.. math::

   D_i^k
   \equiv
   \frac{\sqrt{Pr}}{\sqrt{Ra}} \dder{}{x_j} \dder{u_i^k}{x_j}.

Here the spatial derivatives are approximated by proper finite-difference schemes discussed in :ref:`the other part <spatial_discretisation>`.

**************
Discretisation
**************

A naive temporal discretisation of the momentum equation would be

.. math::

   &\text{do}\,\,\,\, k = 0, 2 \\
   &\,\,\,\,
   \Delta u_i
   =
   \gamma^k \Delta t P_i^{k+1}
   +
   \alpha^k \Delta t g_i^k
   +
   \beta^k  \Delta t g_i^{k-1}, \\
   &\,\,\,\,
   u_i^{k+1}
   =
   u_i^k
   +
   \Delta u_i, \\
   &\text{enddo}

when all advective and diffusive terms are treated explicitly, while

.. math::

   &\text{do}\,\,\,\, k = 0, 2 \\
   &\,\,\,\,
   \Delta u_i
   =
   \gamma^k \Delta t P_i^{k+1}
   +
   \alpha^k \Delta t g_i^k
   +
   \beta^k  \Delta t g_i^{k-1}
   +
   \gamma^k \Delta t h_i^{k  }, \\
   &\,\,\,\,
   u_i^{k+1}
   =
   u_i^k
   +
   \left(
      1 - \frac{\gamma^k \Delta t}{2} \frac{\sqrt{Pr}}{\sqrt{Ra}} \dder{}{y} \dder{}{y}
   \right)^{-1}
   \left(
      1 - \frac{\gamma^k \Delta t}{2} \frac{\sqrt{Pr}}{\sqrt{Ra}} \dder{}{x} \dder{}{x}
   \right)^{-1}
   \Delta u_i, \\
   &\text{enddo}

when diffusive terms are treated implicitly.

Here :math:`g_i^k` and :math:`h_i^k` are used to denote terms which are treated explicitly and implicitly, respectively.
Since the advective terms :math:`A_i^k` are always treated explicitly in time in this project, they are included in :math:`g_i^k`.
Diffusive terms can be fully or partially in :math:`g_i^k` or :math:`h_i^k`, depending on the user specification.

Recall that I regard the three Runge-Kutta sub-steps as the combination of the three Euler-forward schemes having smaller time steps :math:`\gamma^k \Delta t` (see :ref:`the time-marching schemes <time_marchers>`).
As discussed in :ref:`the implicit treatment <implicit_treatment>`, as long as I assume the fluid is incompressible, I need to handle the pressure term implicitly, which is the reason why :math:`P_i^{k+1}` is used.

===========
SMAC method
===========

-------------------------------------------
When diffusive terms are explicitly treated
-------------------------------------------

First, I consider :math:`h_i^k = 0`, i.e. all diffusive terms are treated explicitly.
Compared to the temporal integration of :ref:`the temperature field <temperature_integration>`, I have two additional problems to integrate :ref:`the momentum equation <eq_momentum>`:

   * enforcing the mass conservation while integrating the momentum equation,

   * coupling the pressure, which is independent of time.

One solution to resolve these two challenges is known as Simplified Marker And Cell (SMAC) method (|AMSDEN1970|), which splits the momentum equation into two parts.
In the first step (prediction step), momentum equation is integrated in time, without taking care of the mass conservation.
This velocity field, which in general violates the incompressibility, is corrected by adjusting the pressure field in the second step (correction step).

By introducing a prediction velocity :math:`u_i^*`, this idea to update the velocity field can be written as

.. math::

   &\text{do}\,\,\,\, k = 0, 2 \\
   &\,\,\,\,
   \Delta u_i
   =
   \gamma^k \Delta t P_i^k
   +
   \alpha^k \Delta t g_i^k
   +
   \beta^k  \Delta t g_i^{k-1}, \\
   &\,\,\,\,
   u_i^*
   =
   u_i^k
   +
   \Delta u_i, \\
   &\,\,\,\,
   u_i^{k+1}
   =
   u_i^*
   -
   \gamma^k \Delta t \dder{\psi}{x_i}. \\
   &\text{enddo}

In the above equation, I have :math:`\psi`, which is an unknown scalar potential to be given now.

By adding the above three equations, I have

.. math::

   u_i^{k+1}
   =
   u_i^k
   -
   \gamma^k \Delta t \dder{\left( p^k + \psi \right)}{x_i}
   +
   \alpha^k \Delta t g_i^k
   +
   \beta^k  \Delta t g_i^{k-1}.

As discussed above, since the pressure field is treated implicitly, I obtain

.. math::

   p^{k+1}
   \equiv
   p^{k  }
   +
   \psi.

To compute :math:`\psi` from the known information, I take the (discrete) divergence of the correction step:

.. math::

   u_i^{k+1}
   =
   u_i^*
   -
   \gamma^k \Delta t \dder{\psi}{x_i}

to have

.. math::

   \dder{u_i^{k+1}}{x_i}
   =
   \dder{u_i^{*  }}{x_i}
   -
   \gamma^k \Delta t \dder{}{x_i} \dder{\psi}{x_i}.

By requesting the (discrete) incompressibility constraint on the new velocity field:

.. math::

   \dder{u_i^{k+1}}{x_i}
   =
   0,

I find

.. math::

   \dder{}{x_i} \dder{\psi}{x_i}
   =
   \frac{1}{\gamma^k \Delta t} \dder{u_i^*}{x_i},

which is a Poisson equation with respect to :math:`\psi`.

In summary, the conclusive scheme is

.. math::

   &\text{do}\,\,\,\, k = 0, 2 \\
   &\,\,\,\,
   \Delta u_i
   =
   \gamma^k \Delta t P_i^k
   +
   \alpha^k \Delta t g_i^k
   +
   \beta^k  \Delta t g_i^{k-1}, \\
   &\,\,\,\,
   u_i^*
   =
   u_i^k
   +
   \Delta u_i, \\
   &\,\,\,\,
   \dder{}{x_i} \dder{\psi}{x_i}
   =
   \frac{1}{\gamma^k \Delta t} \dder{u_i^*}{x_i}, \\
   &\,\,\,\,
   u_i^{k+1}
   =
   u_i^*
   -
   \gamma^k \Delta t \dder{\psi}{x_i}, \\
   &\,\,\,\,
   p^{k+1}
   \equiv
   p^{k  }
   +
   \psi. \\
   &\text{enddo}

.. seealso::

   :ref:`fluid/compute_potential <fluid_compute_potential>` for the Poisson solver.

-------------------------------------------
When diffusive terms are implicitly treated
-------------------------------------------

I consider :math:`h_i^k \ne 0`, i.e. all or some diffusive terms are treated implicitly, when minor correction is needed for :math:`\psi`.
The prediction and the correction steps are

.. math::

   &\text{do}\,\,\,\, k = 0, 2 \\
   &\,\,\,\,
   \Delta u_i
   =
   \gamma^k \Delta t P_i^k
   +
   \alpha^k \Delta t g_i^k
   +
   \beta^k  \Delta t g_i^{k-1}
   +
   \gamma^k \Delta t \frac{
      h_i^{k+1}
      +
      h_i^{k  }
   }{2}, \\
   &\,\,\,\,
   u_i^*
   =
   u_i^k
   +
   \Delta u_i, \\
   &\,\,\,\,
   u_i^{k+1}
   =
   u_i^*
   -
   \gamma^k \Delta t \dder{\psi}{x_i}. \\
   &\text{enddo}

From the correction step, I obtain the same Poisson equation

.. math::

   \dder{}{x_i} \dder{\psi}{x_i}
   =
   \frac{1}{\gamma^k \Delta t} \dder{u_i^*}{x_i}.

Since the diffusive terms are treated implicitly now, I have

.. math::

   h_i^k
   =
   D_i^k
   =
   \frac{\sqrt{Pr}}{\sqrt{Ra}} \dder{}{x_j} \dder{u_i^k}{x_j}.

Since :math:`h_i^{k+1}`, which is a function of :math:`u_i^{k+1}`, is unknown in the prediction step, I instead use

.. math::

   \Delta u_i
   =
   \gamma^k \Delta t P_i^k
   +
   \alpha^k \Delta t g_i^k
   +
   \beta^k  \Delta t g_i^{k-1}
   +
   \gamma^k \Delta t \frac{
      h_i^{*  }
      +
      h_i^{k  }
   }{2}

as the prediction step.
Note that :math:`h_i^{*  }` is used instead of :math:`h_i^{k+1}`.

Using the relation coming from the correction step

.. math::

   u_i^*
   =
   u_i^{k+1}
   +
   \gamma^k \Delta t \dder{\psi}{x_i},

I consider to absorb this change in :math:`\psi` by eliminating :math:`u_i^*`:

.. math::

   u_i^{k+1}
   & =
   u_i^k
   -
   \gamma^k \Delta t \dder{}{x_i} \left( \psi + p^k \right)
   +
   \frac{\gamma^k \Delta t^2}{2} \frac{\sqrt{Pr}}{\sqrt{Ra}} \dder{}{x_j} \dder{}{x_j} \gamma^k \dder{\psi}{x_i}
   +
   \cdots \\
   & =
   u_i^k
   -
   \gamma^k \Delta t \dder{p^{k+1}}{x_i}
   +
   \cdots,

giving

.. math::

   \dder{p^{k+1}}{x_i}
   =
   \dder{p^{k  }}{x_i}
   +
   \dder{\psi}{x_i}
   -
   \frac{\gamma^k}{2} \frac{\sqrt{Pr}}{\sqrt{Ra}} \dder{}{x_j} \dder{}{x_j} \Delta t \dder{\psi}{x_i}.

When the discrete gradient operator

.. math::

   \dder{}{x_i}

and the discrete Laplace operator

.. math::

   \dder{}{x_j} \dder{}{x_j}

are commutative, which holds in this project, I have

.. math::

   \dder{p^{k+1}}{x_i}
   =
   \dder{p^{k  }}{x_i}
   +
   \dder{\psi}{x_i}
   -
   \dder{}{x_i} \frac{\gamma^k \Delta t}{2} \frac{\sqrt{Pr}}{\sqrt{Ra}} \dder{}{x_j} \dder{\psi}{x_j}

and thus

.. math::

   p^{k+1}
   =
   p^{k  }
   +
   \psi
   -
   \frac{\gamma^k \Delta t}{2} \frac{\sqrt{Pr}}{\sqrt{Ra}} \dder{}{x_j} \dder{\psi}{x_j}.

In summary, the conclusive scheme is

.. math::

   &\text{do}\,\,\,\, k = 0, 2 \\
   &\,\,\,\,
   \Delta u_i
   =
   \gamma^k \Delta t P_i^k
   +
   \alpha^k \Delta t g_i^k
   +
   \beta^k  \Delta t g_i^{k-1}
   +
   \gamma^k \Delta t h_i^{k  }, \\
   &\,\,\,\,
   u_i^*
   =
   u_i^k
   +
   \left(
      1
      -
      \frac{\gamma^k \Delta t}{2} \frac{\sqrt{Pr}}{\sqrt{Ra}}
      \dder{}{y} \dder{}{y}
   \right)^{-1}
   \left(
      1
      -
      \frac{\gamma^k \Delta t}{2} \frac{\sqrt{Pr}}{\sqrt{Ra}}
      \dder{}{x} \dder{}{x}
   \right)^{-1}
   \Delta u_i, \\
   &\,\,\,\,
   \dder{}{x_i} \dder{\psi}{x_i}
   =
   \frac{1}{\gamma^k \Delta t} \dder{u_i^*}{x_i}, \\
   &\,\,\,\,
   u_i^{k+1}
   =
   u_i^*
   -
   \gamma^k \dder{\psi}{x_i}, \\
   &\,\,\,\,
   p^{k+1}
   \equiv
   p^{k  }
   +
   \psi
   -
   \frac{\gamma^k \Delta t}{2} \frac{\sqrt{Pr}}{\sqrt{Ra}} \dder{}{x_j} \dder{\psi}{x_j}. \\
   &\text{enddo}

.. seealso::

   :ref:`Time-marching schemes <time_marchers>`.

**************
Implementation
**************

#. Compute right-hand-side terms

   :ref:`fluid_compute_rhs <fluid_compute_rhs>`

#. Compute :math:`\Delta u_i`

   :ref:`fluid_predict_field <fluid_predict_field>`

#. Solve linear systems and update velocity field

   :ref:`fluid_predict_field <fluid_predict_field>`

#. Solve Poisson equation

   :ref:`fluid_compute_potential <fluid_compute_potential>`

#. Correct velocity field

   :ref:`fluid_correct_velocity <fluid_correct_velocity>`

#. Update pressure field

   :ref:`fluid_update_pressure.c <fluid_update_pressure>`

