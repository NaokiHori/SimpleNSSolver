
.. _derivation_momentum_diffusive_terms:

###############
Diffusive terms
###############

Diffusive terms in the momentum equations diffuse the momentum in all directions.
They also play a crucial role in diffusing and dissipating the kinetic energy.
In this project, the injected energy should be eventually dissipated by these terms since there is no other energy sink.

******************
Momentum diffusion
******************

Since the diffusive terms are linear (and thus I can treat them implicitly in time, see :ref:`the temporal discretisation <temporal_discretisation>`), I can simply discretise them as

.. math::

   \frac{\sqrt{Pr}}{\sqrt{Ra}} \left\{
      \dder{}{x} \left( \dder{\ux}{x} \right)
      +
      \dder{}{y} \left( \dder{\ux}{y} \right)
   \right\}

in the :math:`x` direction at :math:`\left( \xic, \xjc \right)`, while

.. math::

   \frac{\sqrt{Pr}}{\sqrt{Ra}} \left\{
      \dder{}{x} \left( \dder{\uy}{x} \right)
      +
      \dder{}{y} \left( \dder{\uy}{y} \right)
   \right\}

in the :math:`y` direction at :math:`\left( \yic, \yjc \right)`.
All necessary first-order differentiations are defined, so no interpolation is required.

.. note::

   The results are **not** the same as the corresponding Taylor series expansions of the second-order derivatives when the grid sizes are not equal.

******************
Energy dissipation
******************

In :ref:`the governing equations <governing_equations>`, I took the inner product of :math:`u_i` and the diffusive terms of the momentum equations to derive the relation of the diffusion and the dissipation of :math:`k`:

.. math::

   u_i \der{}{x_j} \left( \der{u_i}{x_j} \right)
   =
   \der{}{x_j} \left( u_i \der{u_i}{x_j} \right)
   -
   \der{u_i}{x_j} \der{u_i}{x_j},

where the pre-factor (the diffusivity) :math:`\sqrt{Pr} / \sqrt{Ra}` are dropped for convenience.
Here I consider the discrete counterpart.

===================
Left-hand-side term
===================

The left-hand-side term yields

.. math::

   \ux
   \left\{
      \dder{}{x} \left( \dder{\ux}{x} \right)
      +
      \dder{}{y} \left( \dder{\ux}{y} \right)
   \right\}
   & =
   \frac{
      \vat{\ux}{\xic, \xjc}
      \vat{
         \dder{\ux}{x}
      }{\xip, \xjc}
      -
      \vat{\ux}{\xic, \xjc}
      \vat{
         \dder{\ux}{x}
      }{\xim, \xjc}
   }{\Delta x_{\xic}} \\
   & +
   \frac{
      \vat{\ux}{\xic, \xjc}
      \vat{
         \dder{\ux}{y}
      }{\xic, \xjp}
      -
      \vat{\ux}{\xic, \xjc}
      \vat{
         \dder{\ux}{y}
      }{\xic, \xjm}
   }{\Delta y}

and

.. math::

   \uy
   \left\{
      \dder{}{x} \left( \dder{\uy}{x} \right)
      +
      \dder{}{y} \left( \dder{\uy}{y} \right)
   \right\}
   & =
   \frac{
      \vat{\uy}{\yic, \yjc}
      \vat{
         \dder{\uy}{x}
      }{\yip, \yjc}
      -
      \vat{\uy}{\yic, \yjc}
      \vat{
         \dder{\uy}{x}
      }{\yim, \yjc}
   }{\Delta x_{\yic}} \\
   & +
   \frac{
      \vat{\uy}{\yic, \yjc}
      \vat{
         \dder{\uy}{y}
      }{\yic, \yjp}
      -
      \vat{\uy}{\yic, \yjc}
      \vat{
         \dder{\uy}{y}
      }{\yic, \yjm}
   }{\Delta y}

at :math:`\left( \xic, \xjc \right)` and :math:`\left( \yic, \yjc \right)`, respectively.

=====================================
The first term in the right-hand side
=====================================

The first term on the right-hand side leads to

.. math::

   \dder{}{x} \left( \ux \dder{\ux}{x} \right)
   +
   \dder{}{y} \left( \ux \dder{\ux}{y} \right)
   & =
   \frac{
      \vat{
         \dintrpa{\ux}{x}
      }{\xip, \xjc}
      \vat{
         \dder{\ux}{x}
      }{\xip, \xjc}
      -
      \vat{
         \dintrpa{\ux}{x}
      }{\xim, \xjc}
      \vat{
         \dder{\ux}{x}
      }{\xim, \xjc}
   }{\Delta x_{\xic}} \\
   & +
   \frac{
      \vat{
         \dintrpa{\ux}{y}
      }{\xic, \xjp}
      \vat{
         \dder{\ux}{y}
      }{\xic, \xjp}
      -
      \vat{
         \dintrpa{\ux}{y}
      }{\xic, \xjm}
      \vat{
         \dder{\ux}{y}
      }{\xic, \xjm}
   }{\Delta y}

and

.. math::

   \dder{}{x} \left( \uy \dder{\uy}{x} \right)
   +
   \dder{}{y} \left( \uy \dder{\uy}{y} \right)
   & =
   \frac{
      \vat{
         \dintrpa{\uy}{x}
      }{\yip, \yjc}
      \vat{
         \dder{\uy}{x}
      }{\yip, \yjc}
      -
      \vat{
         \dintrpa{\uy}{x}
      }{\yim, \yjc}
      \vat{
         \dder{\uy}{x}
      }{\yim, \yjc}
   }{\Delta x_{\yic}} \\
   & +
   \frac{
      \vat{
         \dintrpa{\uy}{y}
      }{\yic, \yjp}
      \vat{
         \dder{\uy}{y}
      }{\yic, \yjp}
      -
      \vat{
         \dintrpa{\uy}{y}
      }{\yic, \yjm}
      \vat{
         \dder{\uy}{y}
      }{\yic, \yjm}
   }{\Delta y}.

.. note::

   For the ease of the later discussion, I keep the above form.
   However, I can further manipulate them to confirm that they are equal to

   .. math::

      \dder{}{x} \dder{}{x} \left( \frac{1}{2} \ux^2 \right)
      +
      \dder{}{y} \dder{}{y} \left( \frac{1}{2} \ux^2 \right)

   and

   .. math::

      \dder{}{x} \dder{}{x} \left( \frac{1}{2} \uy^2 \right)
      +
      \dder{}{y} \dder{}{y} \left( \frac{1}{2} \uy^2 \right),

   which are the discrete diffusive terms of the squared velocity :math:`\ux^2 / 2` and :math:`\uy^2 / 2`, respectively.
   It is readily apparent that they are discretely conservative and thus do not contribute to the net change of the squared velocities, which is consistent with the continuous counterpart.

================
Dissipative term
================

Subtracting the first equation from the second one yields

.. math::

   &
   \frac{
      \left(
         \vat{
            \dintrpa{\ux}{x}
         }{\xip, \xjc}
         -
         \vat{\ux}{\xic, \xjc}
      \right)
      \vat{
         \dder{\ux}{x}
      }{\xip, \xjc}
      -
      \left(
         \vat{
            \dintrpa{\ux}{x}
         }{\xim, \xjc}
         -
         \vat{\ux}{\xic, \xjc}
      \right)
      \vat{
         \dder{\ux}{x}
      }{\xim, \xjc}
   }{\Delta x_{\xic}} \\
   + &
   \frac{
      \left(
         \vat{
            \dintrpa{\ux}{y}
         }{\xic, \xjp}
         -
         \vat{\ux}{\xic, \xjc}
      \right)
      \vat{
         \dder{\ux}{y}
      }{\xic, \xjp}
      -
      \left(
         \vat{
            \dintrpa{\ux}{y}
         }{\xic, \xjm}
         -
         \vat{\ux}{\xic, \xjc}
      \right)
      \vat{
         \dder{\ux}{y}
      }{\xic, \xjm}
   }{\Delta y}

and

.. math::

   &
   \frac{
      \left(
         \vat{
            \dintrpa{\uy}{x}
         }{\yip, \yjc}
         -
         \vat{\uy}{\yic, \yjc}
      \right)
      \vat{
         \dder{\uy}{x}
      }{\yip, \yjc}
      -
      \left(
         \vat{
            \dintrpa{\uy}{x}
         }{\yim, \yjc}
         -
         \vat{\uy}{\yic, \yjc}
      \right)
      \vat{
         \dder{\uy}{x}
      }{\yim, \yjc}
   }{\Delta x_{\yic}} \\
   + &
   \frac{
      \left(
         \vat{
            \dintrpa{\uy}{y}
         }{\yic, \yjp}
         -
         \vat{\uy}{\yic, \yjc}
      \right)
      \vat{
         \dder{\uy}{y}
      }{\yic, \yjp}
      -
      \left(
         \vat{
            \dintrpa{\uy}{y}
         }{\yic, \yjm}
         -
         \vat{\uy}{\yic, \yjc}
      \right)
      \vat{
         \dder{\uy}{y}
      }{\yic, \yjm}
   }{\Delta y}

They are

.. math::

   &
   \frac{1}{2}
   \frac{1}{\Delta x_{\xic}}
   \vat{
      \diffe{\ux}{x}
   }{\xip, \xjc}
   \vat{
      \dder{\ux}{x}
   }{\xip, \xjc}
   +
   \frac{1}{2}
   \frac{1}{\Delta x_{\xic}}
   \vat{
      \diffe{\ux}{x}
   }{\xim, \xjc}
   \vat{
      \dder{\ux}{x}
   }{\xim, \xjc} \\
   + &
   \frac{1}{2}
   \frac{1}{\Delta y}
   \vat{
      \diffe{\ux}{y}
   }{\xic, \xjp}
   \vat{
      \dder{\ux}{y}
   }{\xic, \xjp}
   +
   \frac{1}{2}
   \frac{1}{\Delta y}
   \vat{
      \diffe{\ux}{y}
   }{\xic, \xjm}
   \vat{
      \dder{\ux}{y}
   }{\xic, \xjm} \\
   = &
   \frac{1}{2}
   \frac{\Delta x_{\xip}}{\Delta x_{\xic}}
   \left(
      \vat{
         \dder{\ux}{x}
      }{\xip, \xjc}
   \right)^2
   +
   \frac{1}{2}
   \frac{\Delta x_{\xim}}{\Delta x_{\xic}}
   \left(
      \vat{
         \dder{\ux}{x}
      }{\xim, \xjc}
   \right)^2 \\
   + &
   \frac{1}{2}
   \left(
      \vat{
         \dder{\ux}{y}
      }{\xic, \xjp}
   \right)^2
   +
   \frac{1}{2}
   \left(
      \vat{
         \dder{\ux}{y}
      }{\xic, \xjm}
   \right)^2 \\
   = &
   \color{red}{
      \dintrpv{
         \left( \dder{\ux}{x} \right)^2
      }{x}
      +
      \dintrpa{
         \left( \dder{\ux}{y} \right)^2
      }{y}
   }

and

.. math::

   &
   \frac{1}{2}
   \frac{1}{\Delta x_{\yic}}
   \vat{C}{\yip}
   \vat{
      \diffe{\uy}{x}
   }{\yip, \yjc}
   \vat{
      \dder{\uy}{x}
   }{\yip, \yjc}
   +
   \frac{1}{2}
   \frac{1}{\Delta x_{\yic}}
   \vat{C}{\yim}
   \vat{
      \diffe{\uy}{x}
   }{\yim, \yjc}
   \vat{
      \dder{\uy}{x}
   }{\yim, \yjc} \\
   + &
   \frac{1}{2}
   \frac{1}{\Delta y}
   \vat{
      \diffe{\uy}{y}
   }{\yic, \yjp}
   \vat{
      \dder{\uy}{y}
   }{\yic, \yjp}
   +
   \frac{1}{2}
   \frac{1}{\Delta y}
   \vat{
      \diffe{\uy}{y}
   }{\yic, \yjm}
   \vat{
      \dder{\uy}{y}
   }{\yic, \yjm} \\
   = &
   \frac{1}{2}
   \frac{\Delta x_{\yip}}{\Delta x_{\yic}}
   \vat{C}{\yip}
   \left(
      \vat{
         \dder{\uy}{x}
      }{\yip, \yjc}
   \right)^2
   +
   \frac{1}{2}
   \frac{\Delta x_{\yim}}{\Delta x_{\yic}}
   \vat{C}{\yim}
   \left(
      \vat{
         \dder{\uy}{x}
      }{\yim, \yjc}
   \right)^2 \\
   + &
   \frac{1}{2}
   \left(
      \vat{
         \dder{\uy}{y}
      }{\yic, \yjp}
   \right)^2
   +
   \frac{1}{2}
   \left(
      \vat{
         \dder{\uy}{y}
      }{\yic, \yjm}
   \right)^2 \\
   = &
   \color{red}{
      \dintrpv{
         C \left( \dder{\uy}{x} \right)^2
      }{x}
      +
      \dintrpa{
         \left( \dder{\uy}{y} \right)^2
      }{y}
   }.

Here :math:`C` is a coefficient to correct the boundary values:

.. math::

   C
   =
   \begin{cases}
      \text{wall}      & 2 \\
      \text{otherwise} & 1
   \end{cases}.

This is necessary because :math:`\uy` is defined on the walls and the value is directly used rather than being interpolated.

.. seealso::

   The reddish terms are computed to check :ref:`the instantaneous Nusselt number <nu_kinetic_energy_dissipation_discrete>`.

