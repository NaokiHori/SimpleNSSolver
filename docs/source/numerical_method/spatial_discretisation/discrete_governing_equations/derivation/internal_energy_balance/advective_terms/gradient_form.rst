#############
Gradient form
#############

*******
Summary
*******

I find that the gradient form of the advective terms is

.. math::

   \dintrpv{
      \ux
      \dder{T}{x}
   }{x}
   +
   \dintrpa{
      \uy
      \dder{T}{y}
   }{y}.

**********
Derivation
**********

Analogous to the advective terms in the momentum transfer, the divergence and the gradient forms should be identical to each other via the incompressibility constraint, i.e.

.. math::

   u_i \dder{T}{x_i}
   =
   \dder{u_i T}{x_i}
   -
   T \dder{u_i}{x_i},

which is defined at the cell center :math:`\left( \pic, \pjc \right)`.

Different from the momentumm advection, I do not have to interpolate the incompressibility constraint since it is located at the same position as the temperature.

Thus the gradient form leads

.. math::

   \ux \dder{T}{x}
   & =
   \dder{\ux T}{x}
   -
   T \dder{\ux}{x} \\
   & =
   \dder{
      \ux
      \dintrpa{T}{x}
   }{x}
   - T \dder{
      \ux
   }{x} \\
   & =
   \frac{
      \vat{\ux}{\pip, \pjc}
      \frac{
         \vat{T}{\pipp, \pjc}
         +
         \vat{T}{\pic,  \pjc}
      }{2}
      -
      \vat{\ux}{\pim, \pjc}
      \frac{
         \vat{T}{\pic,  \pjc}
         +
         \vat{T}{\pimm, \pjc}
      }{2}
   }{\Delta x_{\pic}}
   -
   \vat{T}{\pic, \pjc} \frac{
      \vat{\ux}{\pip, \pjc}
      -
      \vat{\ux}{\pim, \pjc}
   }{\Delta x_{\pic}} \\
   & =
   \vat{\ux}{\pip, \pjc}
   \frac{1}{\Delta x_{\pic}}
   \frac{
      \vat{
         \diffe{T}{x}
      }{\pip, \pjc}
   }{2}
   +
   \vat{\ux}{\pim, \pjc}
   \frac{1}{\Delta x_{\pic}}
   \frac{
      \vat{
         \diffe{T}{x}
      }{\pim, \pjc}
   }{2} \\
   & =
   \frac{\Delta x_{\pip}}{2 \Delta x_{\pic}} \vat{
      \left(
         \ux
         \dder{T}{x}
      \right)
   }{\pip, \pjc}
   +
   \frac{\Delta x_{\pim}}{2 \Delta x_{\pic}} \vat{
      \left(
         \ux
         \dder{T}{x}
      \right)
   }{\pim, \pjc} \\
   & =
   \vat{C}{\pip} \vat{
      \left(
         \ux
         \dder{T}{x}
      \right)
   }{\pip, \pjc}
   +
   \vat{C}{\pim} \vat{
      \left(
         \ux
         \dder{T}{x}
      \right)
   }{\pim, \pjc} \\
   & =
   \color{red}{
      \vat{
         \dintrpv{
            \ux
            \dder{T}{x}
         }{x}
      }{\pic, \pjc}
   }

and

.. math::

   \uy \dder{T}{y}
   & =
   \dder{\uy T}{y}
   -
   T \dder{\uy}{y} \\
   & =
   \dder{
      \uy
      \dintrpa{T}{y}
   }{y}
   - T \dder{
      \uy
   }{y} \\
   & =
   \frac{
      \vat{\uy}{\pic, \pjp}
      \frac{
         \vat{T}{\pic, \pjpp}
         +
         \vat{T}{\pic, \pjc}
      }{2}
      -
      \vat{\uy}{\pic, \pjm}
      \frac{
         \vat{T}{\pic, \pjc}
         +
         \vat{T}{\pic, \pjmm}
      }{2}
   }{\Delta y}
   -
   \vat{T}{\pic, \pjc} \frac{
      \vat{\uy}{\pic, \pjp}
      -
      \vat{\uy}{\pic, \pjm}
   }{\Delta y} \\
   & =
   \vat{\uy}{\pic, \pjp}
   \frac{1}{\Delta y}
   \frac{
      \vat{
         \diffe{T}{y}
      }{\pic, \pjp}
   }{2}
   +
   \vat{\uy}{\pic, \pjm}
   \frac{1}{\Delta y}
   \frac{
      \vat{
         \diffe{T}{y}
      }{\pic, \pjm}
   }{2} \\
   & =
   \frac{1}{2} \vat{
      \left(
         \uy
         \dder{T}{y}
      \right)
   }{\pic, \pjp}
   +
   \frac{1}{2} \vat{
      \left(
         \uy
         \dder{T}{y}
      \right)
   }{\pic, \pjm} \\
   & =
   \color{red}{
      \vat{
         \dintrpa{
            \uy
            \dder{T}{y}
         }{y}
      }{\pic, \pjc}
   }.

