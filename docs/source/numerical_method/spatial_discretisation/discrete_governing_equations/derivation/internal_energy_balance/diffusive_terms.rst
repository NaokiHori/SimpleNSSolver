
.. _derivation_internal_diffusive_terms:

###############
Diffusive terms
###############

Diffusive terms in the internal energy equation diffuse the internal energy in all directions.
They also play a crucial role in diffusing and dissipating the quadratic quantity.

********************************
Diffusion of the internal energy
********************************

Since the diffusive terms are linear (and thus I can treat them implicitly in time, see :ref:`the temporal discretisation <temporal_discretisation>`), I can simply discretise them as

.. math::

   \frac{1}{\sqrt{Pr} \sqrt{Ra}} \left\{
      \dder{}{x} \left( \dder{T}{x} \right)
      +
      \dder{}{y} \left( \dder{T}{y} \right)
   \right\}

at :math:`\left( \pic, \pjc \right)` where :math:`T` is defined.
All necessary first-order differentiations are defined, so no interpolation is required.

***************************************************
Diffusion and dissipation of the quadratic quantity
***************************************************

In :ref:`the governing equations <governing_equations>`, I considered the product of :math:`T` and the diffusive terms of the internal energy equation to derive the relation of the diffusion and the dissipation of :math:`h`:

.. math::

   T \der{}{x_i} \left( \der{T}{x_i} \right)
   =
   \der{}{x_i} \left( T \der{T}{x_i} \right)
   -
   \der{T}{x_i} \der{T}{x_i},

where the pre-factor (the diffusivity) :math:`1 / \sqrt{Pr} \sqrt{Ra}` are dropped for convenience.

Here I consider the discrete counterpart.

===================
Left-hand-side term
===================

.. math::

   T
   \left\{
      \dder{}{x} \left( \dder{T}{x} \right)
      +
      \dder{}{y} \left( \dder{T}{y} \right)
   \right\}
   & =
   \frac{
      \vat{T}{\pic, \pjc}
      \vat{
         \dder{T}{x}
      }{\pip, \pjc}
      -
      \vat{T}{\pic, \pjc}
      \vat{
         \dder{T}{x}
      }{\pim, \pjc}
   }{\Delta x_{\pic}} \\
   & +
   \frac{
      \vat{T}{\pic, \pjc}
      \vat{
         \dder{T}{y}
      }{\pic, \pjp}
      -
      \vat{T}{\pic, \pjc}
      \vat{
         \dder{T}{y}
      }{\pic, \pjm}
   }{\Delta y}.

=====================================
The first term on the right-hand side
=====================================

.. math::

   \dder{}{x} \left( T \dder{T}{x} \right)
   +
   \dder{}{y} \left( T \dder{T}{y} \right)
   & =
   \frac{
      \vat{
         \dintrpa{T}{x}
      }{\pip, \pjc}
      \vat{
         \dder{T}{x}
      }{\pip, \pjc}
      -
      \vat{
         \dintrpa{T}{x}
      }{\pim, \pjc}
      \vat{
         \dder{T}{x}
      }{\pim, \pjc}
   }{\Delta x_{\pic}} \\
   & +
   \frac{
      \vat{
         \dintrpa{T}{y}
      }{\pic, \pjp}
      \vat{
         \dder{T}{y}
      }{\pic, \pjp}
      -
      \vat{
         \dintrpa{T}{y}
      }{\pic, \pjm}
      \vat{
         \dder{T}{y}
      }{\pic, \pjm}
   }{\Delta y}.

.. note::

   For the ease of the later discussion, I keep the above form.
   However, I can further manipulate it to confirm that it is equal to

   .. math::

      \dder{}{x} \dder{}{x} \left( \frac{1}{2} T^2 \right)
      +
      \dder{}{y} \dder{}{y} \left( \frac{1}{2} T^2 \right),

   which is the discrete diffusive terms of the squared temperature :math:`T^2 / 2`.
   It is readily apparent that it is discretely conservative and thus does not contribute to the net change of the squared temperature, which is consistent with the continuous counterpart.

================
Dissipative term
================

Subtracting the first equation from the second one yields

.. math::

   &
   \frac{
      \left(
         \vat{
            \dintrpa{T}{x}
         }{\pip, \pjc}
         -
         \vat{T}{\pic, \pjc}
      \right)
      \vat{
         \dder{T}{x}
      }{\pip, \pjc}
      -
      \left(
         \vat{
            \dintrpa{T}{x}
         }{\pim, \pjc}
         -
         \vat{T}{\pic, \pjc}
      \right)
      \vat{
         \dder{T}{x}
      }{\pim, \pjc}
   }{\Delta x_{\pic}} \\
   + &
   \frac{
      \left(
         \vat{
            \dintrpa{T}{y}
         }{\pic, \pjp}
         -
         \vat{T}{\pic, \pjc}
      \right)
      \vat{
         \dder{T}{y}
      }{\pic, \pjp}
      -
      \left(
         \vat{
            \dintrpa{T}{y}
         }{\pic, \pjm}
         -
         \vat{T}{\pic, \pjc}
      \right)
      \vat{
         \dder{T}{y}
      }{\pic, \pjm}
   }{\Delta y} \\
   = &
   \frac{1}{2}
   \frac{1}{\Delta x_{\pic}}
   \vat{C}{\pip}
   \vat{
      \diffe{T}{x}
   }{\pip, \pjc}
   \vat{
      \dder{T}{x}
   }{\pip, \pjc}
   +
   \frac{1}{2}
   \frac{1}{\Delta x_{\pic}}
   \vat{C}{\pim}
   \vat{
      \diffe{T}{x}
   }{\pim, \pjc}
   \vat{
      \dder{T}{x}
   }{\pim, \pjc}
   \\
   + &
   \frac{1}{2}
   \frac{1}{\Delta y}
   \vat{
      \diffe{T}{y}
   }{\pic, \pjp}
   \vat{
      \dder{T}{y}
   }{\pic, \pjp}
   +
   \frac{1}{2}
   \frac{1}{\Delta y}
   \vat{
      \diffe{T}{y}
   }{\pic, \pjm}
   \vat{
      \dder{T}{y}
   }{\pic, \pjm} \\
   = &
   \frac{1}{2}
   \frac{\Delta x_{\pip}}{\Delta x_{\pic}}
   \vat{C}{\pip}
   \left(
      \vat{
         \dder{T}{x}
      }{\pip, \pjc}
   \right)^2
   +
   \frac{1}{2}
   \frac{\Delta x_{\pim}}{\Delta x_{\pic}}
   \vat{C}{\pim}
   \left(
      \vat{
         \dder{T}{x}
      }{\pim, \pjc}
   \right)^2
   \\
   + &
   \frac{1}{2}
   \frac{1}{\Delta y}
   \vat{
      \diffe{T}{y}
   }{\pic, \pjp}
   \vat{
      \dder{T}{y}
   }{\pic, \pjp}
   +
   \frac{1}{2}
   \frac{1}{\Delta y}
   \vat{
      \diffe{T}{y}
   }{\pic, \pjm}
   \vat{
      \dder{T}{y}
   }{\pic, \pjm} \\
   = &
   \color{red}{
      \dintrpv{
         C \left( \dder{T}{x} \right)^2
      }{x}
      +
      \dintrpa{
         \left( \dder{T}{y} \right)^2
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

This is necessary because :math:`T` is defined on the walls and the value is directly used rather than being interpolated.

Now I consider to discretely integrate this equation in the volume.
The dissipative term simply yields

.. math::

   \sum_{\pjc} \sum_{\pic}
   \left[
      \dintrpv{
         C \left( \dder{T}{x} \right)^2
      }{x}
      +
      \dintrpa{
         \left( \dder{T}{y} \right)^2
      }{y}
   \right]
   \Delta x \Delta y,

while the diffusive term yields

.. math::

   \sum_{\pjc} \sum_{\pic}
   \left[
      \dder{}{x} \left( T \dder{T}{x} \right)
      +
      \dder{}{y} \left( T \dder{T}{y} \right)
   \right]
   \Delta x \Delta y
   =
   \sum_{\pjc}
   \left(
      \vat{
         T \dder{T}{x}
      }{x = 1}
      -
      \vat{
         T \dder{T}{x}
      }{x = 0}
   \right)
   \Delta y,

which represents the exchange of the squared temperature through the diffusive process.

.. seealso::

   The reddish terms are computed to check :ref:`the instantaneous Nusselt number <nu_thermal_energy_dissipation_discrete>`.

.. note::

   A more intuitive way to discretise :math:`\der{T}{y}` at :math:`\left( \pic, \pjc \right)` might be

   .. math::

      \frac{
         \vat{T}{\pic, \pjpp}
         -
         \vat{T}{\pic, \pjmm}
      }{2 \Delta y},

   and a similar way in the :math:`x` direction.

   This is :ref:`not conservative <basic_operators_conservative>`, and is clearly different from what is derived in this section.
   As a consequence the energy consistency is broken.
   In particular, the dissipation tends to be underestimated, since the high frequencies are smoothened by the wider stencil.

