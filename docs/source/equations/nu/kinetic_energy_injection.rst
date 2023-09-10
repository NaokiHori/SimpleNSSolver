
.. _nu_kinetic_energy_injection:

####################################################
Nusselt number based on the kinetic energy injection
####################################################

Recall that :ref:`the definition of the Nusselt number <eq_nu_definition>` is

.. math::

   J_{T} \left( x \right)
   =
   Nu
   \times
   J_{T,ref}.

Now I consider to integrate :ref:`the equation of the heat flux <eq_heat_flux>` from the bottom wall :math:`x = 0` to the top wall :math:`x = l_x \equiv 1`, yielding

.. math::

   \int_{x = 0}^{x = 1}
      J_{T} \left( x \right)
   dx
   =
   \int_{x = 0}^{x = 1} Nu J_{T,ref} dx,

or

.. math::

   \int_x \int_z \int_y
      \left[
         \ave{\ux T}{t}
         -
         \frac{1}{\sqrt{Pr} \sqrt{Ra}}
         \ave{\der{T}{x}}{t}
      \right]
   dy dz dx
   =
   Nu J_{T,ref}.

The diffusive terms yield

.. math::

   \int_z \int_y
      \frac{1}{\sqrt{Pr} \sqrt{Ra}}
      \left(
         -
         \vat{\ave{T}{t}}{x = 1}
         +
         \vat{\ave{T}{t}}{x = 0}
      \right)
   dy dz
   =
   \int_z \int_y
      \frac{1}{\sqrt{Pr} \sqrt{Ra}}
   dy dz
   =
   J_{T,ref},

since I impose

.. math::

   \vat{\ave{T}{t}}{x = 0}
   -
   \vat{\ave{T}{t}}{x = 1}
   =
   1

as the boundary conditions.

Thus I notice

.. math::

   \int_x \int_z \int_y
      \ave{\ux T}{t}
   dy dz dx
   +
   J_{T,ref}
   =
   Nu J_{T,ref},

or equivalently

.. math::

   Nu
   =
   \frac{1}{J_{T,ref}}
   \int_V
      \ave{\ux T}{t}
   dV
   +
   1.

Notice that the first term represents the enhancement of the heat transport by the convective effect.
This term appears in :ref:`the equation of the squared velocity <eq_squared_velocity>` as the buoyancy body force, which I interpret as the kinetic energy injection since the squared velocity has a close relation with the kinetic energy.
This is the reason why I call this equation as the Nusselt number based on the kinetic energy injection.

.. note::

   Since I have

   .. math::

      J_{T,ref}
      & =
      \int_z \int_y
         \frac{1}{\sqrt{Pr} \sqrt{Ra}}
      dy dz \\
      & =
      \int_x \int_z \int_y
         \frac{1}{\sqrt{Pr} \sqrt{Ra}}
      dy dz dx \,\, \left( \because l_x \equiv 1 \right) \\
      & =
      \int_V
         \frac{1}{\sqrt{Pr} \sqrt{Ra}}
      dV,

   the result is equal to

   .. math::

      Nu
      & =
      \sqrt{Pr} \sqrt{Ra}
      \frac{
         \int_V
            \ave{\ux T}{t}
         dV
      }{
         \int_V dV
      }
      +
      1 \\
      & =
      \sqrt{Pr} \sqrt{Ra} \ave{\ux T}{x,y,z,t}
      +
      1,

   which might be more general.

.. seealso::

   :ref:`Discrete counterpart <nu_kinetic_energy_injection_discrete>`.

