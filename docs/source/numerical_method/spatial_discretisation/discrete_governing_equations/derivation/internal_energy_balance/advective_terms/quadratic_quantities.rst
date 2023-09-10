####################
Quadratic quantities
####################

*******
Summary
*******

I find the discrete form of the advection of :math:`T` is given by

.. math::

   \dder{\ux \dintrpa{T}{x}}{x}
   +
   \dder{\uy \dintrpa{T}{y}}{y},

where the interpolation in the :math:`x` direction which was undefined is explicitly given.

**********
Derivation
**********

I consider the advection of the quadratic quantity :math:`h`.
As shown in :ref:`the governing equations <governing_equations>`, this is obtained by multiplying the advective terms of the internal energy advection by :math:`T`.
The discrete counterpart is twofold,

.. math::

   T \dder{\ux \dintrpu{T}{x}}{x}
   =
   \color{blue}{
   \vat{
      T
   }{\pic, \pjc}
   }
   \frac{
      \color{blue}{
      \vat{\ux}{\pip, \pjc}
      }
      \left(
         \color{blue}{
         \vat{c^+}{\pipp} \vat{T}{\pipp, \pjc}
         }
         +
         \vat{c^+}{\pic } \vat{T}{\pic,  \pjc}
      \right)
      -
      \color{blue}{
      \vat{\ux}{\pim, \pjc}
      }
      \left(
         \vat{c^-}{\pic } \vat{T}{\pic,  \pjc}
         +
         \color{blue}{
         \vat{c^-}{\pimm} \vat{T}{\pimm, \pjc}
         }
      \right)
   }{\Delta x_{\pic}}

and

.. math::

   T \dder{\uy \dintrpu{T}{y}}{y}
   =
   \vat{
      T
   }{\pic, \pjc}
   \frac{
      \vat{\uy}{\pic, \pjp}
      \frac{
         \vat{T}{\pic, \pjpp}
         +
         \vat{T}{\pic, \pjc }
      }{2}
      -
      \vat{\uy}{\pic, \pjm}
      \frac{
         \vat{T}{\pic, \pjc }
         +
         \vat{T}{\pic, \pjmm}
      }{2}
   }{\Delta y},

where small letters :math:`c` are coefficients to be determined.

To keep the net amount of :math:`h`, I enforce them to be conservative.
The second part yields

.. math::

   \dder{\uy q_y}{y}
   +
   \frac{
      \vat{T^2}{\pic, \pjc}
   }{2}
   \dder{\uy}{y},

where :math:`q_y` is the quadratic quantity, defined as the product of the surrounding :math:`T` at cell faces where :math:`\uy` are located, e.g.

.. math::

   \vat{q_y}{\pip, \pjc}
   \equiv
   \frac{1}{2}
   \vat{T}{\pic,  \pjc}
   \vat{T}{\pipp, \pjc}

at :math:`\left( \pip, \pjc \right)`.

Since the first term is in a conservative form, the volume integral vanishes, while the second part is the residual, which is generally non-zero.

Similarly, I have the following two terms in the :math:`x` component.

==================
Quadratic quantity
==================

To be consistent with the :math:`y` component, the bluish terms should include the other quadratic quantity :math:`q_x` and I expect them to be conservative.
To achieve this, I request the coefficients :math:`\vat{c^+}{\pipp}` and :math:`\vat{c^-}{\pimm}` to be :math:`1/2`.

========
Residual
========

The second part leads to

.. math::

   \vat{T^2}{\pic, \pjc}
   \frac{
      \vat{c^+}{\pic }
      \vat{\ux}{\pim, \pjc}
      -
      \vat{c^-}{\pic }
      \vat{\ux}{\pim, \pjc}
   }{\Delta x_{\pic}}.

This should be canceled out with the :math:`y` residual:

.. math::

   \frac{
      \vat{T^2}{\pic, \pjc}
   }{2}
   \dder{\uy}{y},

concluding that the two coefficients :math:`\vat{c^+}{\pic }` and :math:`\vat{c^-}{\pic }` to be :math:`1/2`.

In summary, all undefined coefficients are :math:`1/2` and thus I find :math:`\dintrpu{T}{x}` is the arithmetic average :math:`\dintrpa{T}{x}`.

Finally I obtain the advective terms in divergence form:

.. math::

   \dder{\ux T}{x}
   +
   \dder{\uy T}{y}
   =
   \color{red}{
      \dder{
         \ux
         \dintrpa{T}{x}
      }{x}
      +
      \dder{
         \uy
         \dintrpa{T}{y}
      }{y}
   }

defined at the cell center :math:`\left( \pic, \pjc \right)`.

