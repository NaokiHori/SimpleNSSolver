#############
Gradient form
#############

*******
Summary
*******

Finally I find that the gradient form of the advective terms can be written as

.. math::

   \dintrpv{
      \dintrpa{\ux}{x}
      \dder{\ux}{x}
   }{x}
   +
   \dintrpa{
      \dintrpv{\uy}{x}
      \dder{\ux}{y}
   }{y}

at :math:`\left( \xic, \xjc \right)` and

.. math::

   \dintrpv{
      \dintrpa{\ux}{y}
      \dder{\uy}{x}
   }{x}
   +
   \dintrpa{
      \dintrpa{\uy}{y}
      \dder{\uy}{y}
   }{y}

at :math:`\left( \yic, \yjc \right)`.

**********
Derivation
**********

The gradient form of the advective terms in the :math:`x` momentum equation

.. math::

   \ux \dder{\ux}{x}
   +
   \uy \dder{\ux}{y},

and the terms in the :math:`y` momentum equation

.. math::

   \ux \dder{\uy}{x}
   +
   \uy \dder{\uy}{y}

are written as

.. math::

   \dder{\ux \ux}{x}
   +
   \dder{\uy \ux}{y}
   -
   \ux \left(
      \dder{\ux}{x}
      +
      \dder{\uy}{y}
   \right)

and

.. math::

   \dder{\ux \uy}{x}
   +
   \dder{\uy \uy}{y}
   -
   \uy \left(
      \dder{\ux}{x}
      +
      \dder{\uy}{y}
   \right),

i.e.

.. math::

   \left( \text{divergence form} \right) - \ux \left( \text{incompressibility constraint} \right),

and

.. math::

   \left( \text{divergence form} \right) - \uy \left( \text{incompressibility constraint} \right).

Since the incompressibility constraint is defined at each cell center, I need to interpolate it to the cell face in-between:

.. math::

   \vat{
      \dder{u_i}{x_i}
   }{\xic, \xjc}
   =
   \vat{
      \dintrpu{
         \dder{u_i}{x_i}
      }{x}
   }{\xic, \xjc},

and

.. math::

   \vat{
      \dder{u_i}{x_i}
   }{\yic, \yjc}
   =
   \vat{
      \dintrpa{
         \dder{u_i}{x_i}
      }{y}
   }{\yic, \yjc}.

In the :math:`y` direction, the interpolation is already known.
In the :math:`x` direction, I need to choose a proper interpolation.

.. note::

   :math:`\vat{\ux}{\xic,\xjc}` in front of the incompressibility constraint is controversial since there are mainly two possibilities:

   #. Use local value at :math:`\left( \xic, \xjc \right)`

      .. math::

         \vat{
            \left(
               \ux
               \dintrpu{
                  \dder{u_i}{x_i}
               }{x}
            \right)
         }{\xic, \xjc}

   #. Interpolate as well as the incompressibility constraint

      .. math::

         \vat{
            \dintrpu{
               \left(
                  \dintrpa{\ux}{x}
                  \dder{u_i}{x_i}
               \right)
            }{x}
         }{\xic, \xjc}

   To answer this question, I focus on the :math:`y` component :math:`\uy \der{\ux}{y}`.
   As shown in the divergence form, three :math:`\ux` are involved: :math:`\vat{\ux}{\xic,\xjp}`, :math:`\vat{\ux}{\xic,\xjc}`, and :math:`\vat{\ux}{\xic,\xjm}`, which are defined on the same :math:`x` cell face.

   If I adopt the second option, two additional :math:`\ux` components located at different :math:`x` cell faces are involved, which is not consistent because the divergence and the gradient forms have different :math:`\ux` information.

   Thus I notice that I should take the first option.

To go further, I let

.. math::

   \vat{\dintrpu{q}{x}}{\xic}
   =
   \vat{C}{\xip} \vat{q}{\xip}
   +
   \vat{C}{\xim} \vat{q}{\xim}.

Since this is a linear operation, the interpolations and the differentiations should be interchangeable:

.. math::

   \dintrpu{
      \dder{\ux}{x}
   }{x}
   =
   \dder{
      \dintrpu{\ux}{x}
   }{x},
   \dintrpu{
      \dder{\uy}{y}
   }{x}
   =
   \dder{
      \dintrpu{\uy}{x}
   }{y}.

The second relation is obvious since the operations in the :math:`x` and the :math:`y` directions are independent.

To satisfy the first relation, I should properly determine :math:`\vat{C}{\xip}` and :math:`\vat{C}{\xim}`.
The interpolation of the differentiation leads to

.. math::

   \vat{
      \dintrpu{
         \dder{\ux}{x}
      }{x}
   }{\xic, \xjc}
   & =
   \vat{C}{\xip} \vat{\dder{\ux}{x}}{\xip, \xjc}
   +
   \vat{C}{\xim} \vat{\dder{\ux}{x}}{\xim, \xjc} \\
   & =
   \vat{C}{\xip} \frac{
      \vat{\ux}{\xipp, \xjc}
      -
      \vat{\ux}{\xic,  \xjc}
   }{\Delta x_{\xip}}
   +
   \vat{C}{\xim} \frac{
      \vat{\ux}{\xic,  \xjc}
      -
      \vat{\ux}{\ximm, \xjc}
   }{\Delta x_{\xim}},

while the differentiation of the interpolation leads to

.. math::

   \vat{
      \dder{
         \dintrpu{\ux}{x}
      }{x}
   }{\xic, \xjc}
   & =
   \frac{
      \vat{\dintrpu{\ux}{x}}{\xip, \xjc}
      -
      \vat{\dintrpu{\ux}{x}}{\xim, \xjc}
   }{\Delta x_{\xic}} \\
   & =
   \frac{
      \vat{\dintrpa{\ux}{x}}{\xip, \xjc}
      -
      \vat{\dintrpa{\ux}{x}}{\xim, \xjc}
   }{\Delta x_{\xic}} \\
   & =
   \frac{
      \frac{
         \vat{\ux}{\xipp, \xjc}
         +
         \vat{\ux}{\xic,  \xjc}
      }{2}
      -
      \frac{
         \vat{\ux}{\xic,  \xjc}
         +
         \vat{\ux}{\ximm, \xjc}
      }{2}
   }{\Delta x_{\xic}}.

By comparing these two equations, I notice

.. math::

   \vat{C}{\xip}
   & =
   \frac{\Delta x_{\xip}}{2 \Delta x_{\xic}}, \\
   \vat{C}{\xim}
   & =
   \frac{\Delta x_{\xim}}{2 \Delta x_{\xic}},

which is the volume average: :math:`\dintrpv{q}{x}`.

Thus the gradient form of the advective terms is

.. math::

   \ux \dder{\ux}{x}
   & =
   \dder{\ux \ux}{x}
   -
   \ux \dder{\ux}{x} \\
   & =
   \dder{
      \dintrpa{\ux}{x}
      \dintrpa{\ux}{x}
   }{x}
   -
   \ux \dintrpv{
      \dder{\ux}{x}
   }{x} \\
   & =
   \dder{
      \dintrpa{\ux}{x}
      \dintrpa{\ux}{x}
   }{x}
   -
   \ux \dder{
      \dintrpa{\ux}{x}
   }{x} \\
   & =
   \frac{
      \vat{\left(
         \dintrpa{\ux}{x}
         \dintrpa{\ux}{x}
      \right)}{\xip, \xjc}
      -
      \vat{\left(
         \dintrpa{\ux}{x}
         \dintrpa{\ux}{x}
      \right)}{\xim, \xjc}
   }{\Delta x_{\xic}}
   -
   \vat{\ux}{\xic, \xjc}
   \frac{
      \vat{\dintrpa{\ux}{x}}{\xip, \xjc}
      -
      \vat{\dintrpa{\ux}{x}}{\xim, \xjc}
   }{\Delta x_{\xic}} \\
   & =
   \vat{\dintrpa{\ux}{x}}{\xip, \xjc}
   \frac{
      \vat{\dintrpa{\ux}{x}}{\xip, \xjc}
      -
      \vat{\ux}{\xic, \xjc}
   }{\Delta x_{\xic}}
   -
   \vat{\dintrpa{\ux}{x}}{\xim, \xjc}
   \frac{
      \vat{\dintrpa{\ux}{x}}{\xim, \xjc}
      -
      \vat{\ux}{\xic, \xjc}
   }{\Delta x_{\xic}} \\
   & =
   \vat{\dintrpa{\ux}{x}}{\xip, \xjc}
   \frac{1}{\Delta x_{\xic}} \frac{
      \vat{\diffe{\ux}{x}}{\xip, \xjc}
   }{2}
   +
   \vat{\dintrpa{\ux}{x}}{\xim, \xjc}
   \frac{1}{\Delta x_{\xic}} \frac{
      \vat{\diffe{\ux}{x}}{\xim, \xjc}
   }{2} \\
   & =
   \frac{\Delta x_{\xip}}{2 \Delta x_{\xic}}
   \vat{\left( \dintrpa{\ux}{x} \dder{\ux}{x} \right)}{\xip, \xjc}
   +
   \frac{\Delta x_{\xim}}{2 \Delta x_{\xic}}
   \vat{\left( \dintrpa{\ux}{x} \dder{\ux}{x} \right)}{\xim, \xjc} \\
   & =
   \vat{C}{\xip}
   \vat{\left( \dintrpa{\ux}{x} \dder{\ux}{x} \right)}{\xip, \xjc}
   +
   \vat{C}{\xim}
   \vat{\left( \dintrpa{\ux}{x} \dder{\ux}{x} \right)}{\xim, \xjc} \\
   & =
   \color{red}{
      \vat{
         \dintrpv{
            \dintrpa{\ux}{x} \dder{\ux}{x}
         }{x}
      }{\xic, \xjc}
   },

.. math::

   \uy \dder{\ux}{x}
   & =
   \dder{\uy \ux}{x}
   -
   \uy \dder{\ux}{x} \\
   & =
   \dder{
      \dintrpv{\uy}{x}
      \dintrpa{\ux}{y}
   }{y}
   -
   \ux \dintrpv{
      \dder{\uy}{y}
   }{x} \\
   & =
   \dder{
      \dintrpv{\uy}{x}
      \dintrpa{\ux}{y}
   }{y}
   -
   \ux \dder{
      \dintrpv{\uy}{x}
   }{y} \\
   & =
   \frac{
      \vat{
         \left(
            \dintrpv{\uy}{x}
            \dintrpa{\ux}{y}
         \right)
      }{\xic, \xjp}
      -
      \vat{
         \left(
            \dintrpv{\uy}{x}
            \dintrpa{\ux}{y}
         \right)
      }{\xic, \xjm}
   }{\Delta y}
   -
   \vat{
      \ux
   }{\xic, \xjc}
   \frac{
      \vat{
         \dintrpv{\uy}{x}
      }{\xic, \xjp}
      -
      \vat{
         \dintrpv{\uy}{x}
      }{\xic, \xjm}
   }{\Delta y} \\
   & =
   \vat{\dintrpv{\uy}{x}}{\xic, \xjp}
   \frac{
      \vat{\dintrpa{\ux}{y}}{\xic, \xjp}
      -
      \vat{\ux}{\xic, \xjc}
   }{\Delta y}
   -
   \vat{\dintrpv{\uy}{x}}{\xic, \xjm}
   \frac{
      \vat{\dintrpa{\ux}{y}}{\xic, \xjm}
      -
      \vat{\ux}{\xic, \xjc}
   }{\Delta y} \\
   & =
   \frac{1}{2}
   \vat{
      \left(
         \dintrpv{\uy}{x}
         \dder{\ux}{y}
      \right)
   }{\xic, \xjp}
   +
   \frac{1}{2} \vat{
      \left(
         \dintrpv{\uy}{x}
         \dder{\ux}{y}
      \right)
   }{\xic, \xjm} \\
   & =
   \color{red}{
      \vat{
         \dintrpa{
            \dintrpv{\uy}{x}
            \dder{\ux}{y}
         }{y}
      }{\xic, \xjc}
   },

.. math::

   \ux \dder{\uy}{x}
   & =
   \dder{\ux \uy}{x}
   -
   \uy \dder{\ux}{x} \\
   & =
   \dder{
      \dintrpa{\ux}{y}
      \dintrpa{\uy}{x}
   }{x}
   -
   \uy \dintrpa{
      \dder{\ux}{x}
   }{y} \\
   & =
   \dder{
      \dintrpa{\ux}{y}
      \dintrpa{\uy}{x}
   }{x}
   -
   \uy \dder{
      \dintrpa{\ux}{y}
   }{x} \\
   & =
   \frac{
      \vat{\left(
         \dintrpa{\ux}{y}
         \dintrpa{\uy}{x}
      \right)}{\yip, \yjc}
      -
      \vat{\left(
         \dintrpa{\ux}{y}
         \dintrpa{\uy}{x}
      \right)}{\yim, \yjc}
   }{\Delta x_{\yic}}
   -
   \vat{\uy}{\yic, \yjc}
   \frac{
      \vat{\dintrpa{\ux}{y}}{\yip, \yjc}
      -
      \vat{\dintrpa{\ux}{y}}{\yim, \yjc}
   }{\Delta x_{\yic}} \\
   & =
   \vat{\dintrpa{\ux}{y}}{\yip, \yjc}
   \frac{
      \vat{\dintrpa{\uy}{x}}{\yip, \yjc}
      -
      \vat{\uy}{\yic, \yjc}
   }{\Delta x_{\yic}}
   -
   \vat{\dintrpa{\ux}{y}}{\yim, \yjc}
   \frac{
      \vat{\dintrpa{\uy}{x}}{\yim, \yjc}
      -
      \vat{\uy}{\yic, \yjc}
   }{\Delta x_{\yic}} \\
   & =
   \vat{\dintrpa{\ux}{y}}{\yip, \yjc}
   \frac{1}{\Delta x_{\yic}} \frac{
      \vat{\diffe{\uy}{x}}{\yip, \yjc}
   }{2}
   +
   \vat{\dintrpa{\ux}{y}}{\yim, \yjc}
   \frac{1}{\Delta x_{\yic}} \frac{
      \vat{\diffe{\uy}{x}}{\yim, \yjc}
   }{2} \\
   & =
   \frac{\Delta x_{\yip}}{2 \Delta x_{\yic}}
   \vat{\left( \dintrpa{\ux}{y} \dder{\uy}{x} \right)}{\yip, \yjc}
   +
   \frac{\Delta x_{\yim}}{2 \Delta x_{\yic}}
   \vat{\left( \dintrpa{\ux}{y} \dder{\uy}{x} \right)}{\yim, \yjc} \\
   & =
   \vat{C}{\yip}
   \vat{\left( \dintrpa{\ux}{y} \dder{\uy}{x} \right)}{\yip, \yjc}
   +
   \vat{C}{\yim}
   \vat{\left( \dintrpa{\ux}{y} \dder{\uy}{x} \right)}{\yim, \yjc} \\
   & =
   \color{red}{
      \vat{
         \dintrpv{
            \dintrpa{\ux}{y} \dder{\uy}{x}
         }{x}
      }{\yic, \yjc}
   },

and

.. math::
   \uy \dder{\uy}{y}
   & =
   \dder{\uy \uy}{y}
   -
   \uy \dder{\uy}{y} \\
   & =
   \dder{
      \dintrpa{\uy}{y}
      \dintrpa{\uy}{y}
   }{y}
   -
   \uy \dintrpa{
      \dder{\uy}{y}
   }{y} \\
   & =
   \dder{
      \dintrpa{\uy}{y}
      \dintrpa{\uy}{y}
   }{y}
   -
   \uy \dder{
      \dintrpa{\uy}{y}
   }{y} \\
   & =
   \frac{
      \vat{
         \left(
            \dintrpa{\uy}{y}
            \dintrpa{\uy}{y}
         \right)
      }{\yic, \yjp}
      -
      \vat{
         \left(
            \dintrpa{\uy}{y}
            \dintrpa{\uy}{y}
         \right)
      }{\yic, \yjm}
   }{\Delta y}
   -
   \vat{
      \uy
   }{\yic, \yjc}
   \frac{
      \vat{
         \dintrpa{\uy}{y}
      }{\yic, \yjp}
      -
      \vat{
         \dintrpa{\uy}{y}
      }{\yic, \yjm}
   }{\Delta y} \\
   & =
   \vat{\dintrpa{\uy}{y}}{\yic, \yjp}
   \frac{
      \vat{\dintrpa{\uy}{y}}{\yic, \yjp}
      -
      \vat{\uy}{\yic, \yjc}
   }{\Delta y}
   -
   \vat{\dintrpa{\uy}{y}}{\yic, \yjm}
   \frac{
      \vat{\dintrpa{\uy}{y}}{\yic, \yjm}
      -
      \vat{\uy}{\yic, \yjc}
   }{\Delta y} \\
   & =
   \frac{1}{2} \vat{
      \left(
         \dintrpa{\uy}{y}
         \dder{\uy}{y}
      \right)
   }{\yic, \yjp}
   +
   \frac{1}{2} \vat{
      \left(
         \dintrpa{\uy}{y}
         \dder{\uy}{y}
      \right)
   }{\yic, \yjm} \\
   & =
   \color{red}{
      \vat{
         \dintrpa{
            \dintrpa{\uy}{y}
            \dder{\uy}{y}
         }{y}
      }{\yic, \yjc}
   }.

*******************
Energy conservation
*******************

Since the gradient form is equivalent to the divergence form as long as the incompressibility constaint is satisfied, the gradient form should conserve the discrete kinetic energy.
Here, I confirm this fact just for completeness.
By multiplying the :math:`x` component

.. math::

   \vat{
      \dintrpv{
         \dintrpa{\ux}{x} \dder{\ux}{x}
      }{x}
   }{\xic, \xjc}
   +
   \vat{
      \dintrpa{
         \dintrpv{\uy}{x}
         \dder{\ux}{y}
      }{y}
   }{\xic, \xjc}

with :math:`\vat{\ux}{\xic, \xjc}`, I obtain

.. math::

   & \vat{\ux}{\xic, \xjc}
   \left(
      \vat{
         \dintrpv{
            \dintrpa{\ux}{x} \dder{\ux}{x}
         }{x}
      }{\xic, \xjc}
      +
      \vat{
         \dintrpa{
            \dintrpv{\uy}{x}
            \dder{\ux}{y}
         }{y}
      }{\xic, \xjc}
   \right) \\
   & =
   \vat{\ux}{\xic, \xjc}
   \frac{\Delta x_{\xip}}{2 \Delta x_{\xic}}
   \vat{\left( \dintrpa{\ux}{x} \dder{\ux}{x} \right)}{\xip, \xjc}
   +
   \vat{\ux}{\xic, \xjc}
   \frac{\Delta x_{\xim}}{2 \Delta x_{\xic}}
   \vat{\left( \dintrpa{\ux}{x} \dder{\ux}{x} \right)}{\xim, \xjc} \\
   & +
   \vat{\ux}{\xic, \xjc}
   \frac{1}{2}
   \vat{
      \left(
         \dintrpv{\uy}{x}
         \dder{\ux}{y}
      \right)
   }{\xic, \xjp}
   +
   \vat{\ux}{\xic, \xjc}
   \frac{1}{2} \vat{
      \left(
         \dintrpv{\uy}{x}
         \dder{\ux}{y}
      \right)
   }{\xic, \xjm}.

By rearranging the differentiations (while keeping the averages), I have

.. math::

   & \dder{
      \dintrpa{\ux}{x}
      q_{xx}
   }{x}
   -
   \frac{
      \vat{\ux^2}{\xic, \xjc}
   }{2}
   \frac{
      \vat{
         \dintrpa{\ux}{x}
      }{\xip, \xjc}
      -
      \vat{
         \dintrpa{\ux}{x}
      }{\xim, \xjc}
   }{\Delta x_{\xic}} \\
   + &
   \dder{
      \dintrpu{\uy}{x}
      q_{xy}
   }{y}
   -
   \frac{
      \vat{\ux^2}{\xic, \xjc}
   }{2}
   \frac{
      \vat{
         \dintrpu{\uy}{x}
      }{\xic, \xjp}
      -
      \vat{
         \dintrpu{\uy}{x}
      }{\xic, \xjm}
   }{\Delta y},

which are consisted of the conservative terms and the volume-weighted continuities at :math:`\left( \xim, \xjc \right)` and :math:`\left( \xip, \xjc \right)` and thus vanish when integrated in the whole domain.

Similarly, the advective terms in the :math:`y` direction in the gradient form yields

.. math::

   & \vat{\uy}{\yic, \yjc}
   \left(
      \vat{
         \dintrpv{
            \dintrpa{\ux}{y} \dder{\uy}{x}
         }{x}
      }{\yic, \yjc}
      +
      \vat{
         \dintrpa{
            \dintrpa{\uy}{y}
            \dder{\uy}{y}
         }{y}
      }{\yic, \yjc}
   \right) \\
   & =
   \vat{\uy}{\yic, \yjc}
   \frac{\Delta x_{\yip}}{2 \Delta x_{\yic}}
   \vat{\left( \dintrpa{\ux}{y} \dder{\uy}{x} \right)}{\yip, \yjc}
   +
   \vat{\uy}{\yic, \yjc}
   \frac{\Delta x_{\yim}}{2 \Delta x_{\yic}}
   \vat{\left( \dintrpa{\ux}{y} \dder{\uy}{x} \right)}{\yim, \yjc} \\
   & +
   \vat{\uy}{\yic, \yjc}
   \frac{1}{2} \vat{
      \left(
         \dintrpa{\uy}{y}
         \dder{\uy}{y}
      \right)
   }{\yic, \yjp}
   +
   \vat{\uy}{\yic, \yjc}
   \frac{1}{2} \vat{
      \left(
         \dintrpa{\uy}{y}
         \dder{\uy}{y}
      \right)
   }{\yic, \yjm}.

Again, I have the conservative terms and the volume-averaged continuities at :math:`\left( \yic, \yjm \right)` and :math:`\left( \yic, \yjp \right)`, which vanish when integrated in the whole domain.

