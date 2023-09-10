##################
Quadratic quantity
##################

*******
Summary
*******

In this part, I show the discrete advective terms of the momentum equations are give by

.. math::

   \dder{
      \dintrpa{\ux}{x}
      \dintrpa{\ux}{x}
   }{x}
   +
   \dder{
      \dintrpv{\uy}{x}
      \dintrpa{\ux}{y}
   }{y}

at :math:`\left( \xic, \xjc \right)` and

.. math::

   \dder{
      \dintrpa{\ux}{y}
      \dintrpa{\uy}{x}
   }{x}
   +
   \dder{
      \dintrpa{\uy}{y}
      \dintrpa{\uy}{y}
   }{y}

at :math:`\left( \yic, \yjc \right)`, respectively.

Note that the interpolations :math:`\dintrpv{\uy}{x}` and :math:`\dintrpa{\ux}{y}`, which were undefined in :ref:`the previous part <momentum_advective_terms_divergence_form>`, are now concluded.

Also, the discrete analogue of the advection of the kinetic energy leads to

.. math::

   \dder{
      \dintrpa{\ux}{x}
      q_{xx}
   }{x}
   +
   \dder{
      \dintrpv{\uy}{x}
      q_{xy}
   }{y}

and

.. math::

   \dder{
      \dintrpa{\ux}{y}
      q_{yx}
   }{x}
   +
   \dder{
      \dintrpa{\uy}{y}
      q_{yy}
   }{y},

where :math:`q_{ij}` are the discrete analogues of the squared velocity defined as

.. math::

   \vat{q_{xx}}{\xip, \xjc}
   \equiv
   \frac{1}{2}
   \vat{\ux}{\xipp, \xjc}
   \vat{\ux}{\xic,  \xjc},

.. math::

   \vat{q_{xy}}{\xic, \xjp}
   \equiv
   \frac{1}{2}
   \vat{\ux}{\xic, \xjpp}
   \vat{\ux}{\xic, \xjc },

and

.. math::

   \vat{q_{yx}}{\yip, \yjc}
   \equiv
   \frac{1}{2}
   \vat{\uy}{\yipp, \yjc}
   \vat{\uy}{\yic,  \yjc},

.. math::

   \vat{q_{yy}}{\yic, \yjp}
   \equiv
   \frac{1}{2}
   \vat{\uy}{\yic, \yjpp}
   \vat{\uy}{\yic, \yjc }.

An important takeaway here is that these terms are conservative and vanish when integrated in the whole volume:

.. math::

   \sum_{\xjc} \sum_{\xic}
   \left[
      \dder{
         \dintrpa{\ux}{x}
         q_{xx}
      }{x}
      +
      \dder{
         \dintrpv{\uy}{x}
         q_{xy}
      }{y}
   \right]_{\xic, \xjc}
   \Delta x_{\xic} \Delta y_{\xjc}
   =
   0,

.. math::

   \sum_{\yjc} \sum_{\yic}
   \left[
      \dder{
         \dintrpa{\ux}{y}
         q_{yx}
      }{x}
      +
      \dder{
         \dintrpa{\uy}{y}
         q_{yy}
      }{y}
   \right]_{\yic, \yjc}
   \Delta x_{\yic} \Delta y_{\yjc}
   =
   0,

i.e. they do not contribute to the increase nor decrease in the squared velocities.

**********
Derivation
**********

In the continuous domain, the advection of the kinetic energy is obtained by taking the inner product of the velocity and the momentum advections.
Here I consider the corresponding relation in the discrete domain.

.. note::

   In the continuous space, the local kinetic energy

   .. math::

      k = k \left( t, x, y \right)

   is clearly the sum of

   .. math::

      \frac{1}{2} \left[ \ux \left( t, x, y \right) \right]^2

   and

   .. math::

      \frac{1}{2} \left[ \uy \left( t, x, y \right) \right]^2.

   After discretised, however, the definition of the kinetic energy becomes ambiguous since :math:`\ux` and :math:`\uy` are defined at different positions.

   In this project, I consider these two contributions separately:

   #. :math:`\ux` contribution

      .. math::

         \der{}{t} \vat{
            \left( \frac{1}{2} \ux^2 \right)
         }{\xic, \xjc} + \cdots,

      which is obtained by multiplying the discrete :math:`x` momentum balance by :math:`\ux`.

   #. :math:`\uy` contribution

      .. math::

         \der{}{t} \vat{
            \left( \frac{1}{2} \uy^2 \right)
         }{\yic, \yjc} + \cdots,

      which is obtained by multiplying the discrete :math:`y` momentum balance by :math:`\uy`.

   I define the total and the discrete kinetic energy :math:`K = K \left( t \right)` as the sum of the discrete volume integrals of these two equations:

   .. math::

      \der{K}{t}
      =
      \sum_{\xic} \sum_{\xjc} \der{}{t} \vat{
         \left(
            \frac{1}{2} \ux^2
            \Delta x \Delta y
         \right)
      }{\xic, \xjc}
      +
      \sum_{\yic} \sum_{\yjc} \der{}{t} \vat{
         \left(
            \frac{1}{2} \uy^2
            \Delta x \Delta y
         \right)
      }{\yic, \yjc},

   which is investigated in this section.

At :math:`\left( \xic, \xjc \right)` where :math:`\ux` is defined, by multiplying the :math:`x` momentum by :math:`\ux`, I have

.. math::

   \ux \left(
      \dder{
         \dintrpa{\ux}{x}
         \dintrpa{\ux}{x}
      }{x}
      +
      \dder{
         \dintrpu{\uy}{x}
         \dintrpa{\ux}{y}
      }{y}
   \right).

The first term leads to

.. math::

   & \vat{\ux}{\xic, \xjc}
   \frac{
      \vat{
         \dintrpa{\ux}{x}
      }{\xip, \xjc}
      \frac{
         \vat{\ux}{\xipp, \xjc}
         +
         \vat{\ux}{\xic , \xjc}
      }{2}
      -
      \vat{
         \dintrpa{\ux}{x}
      }{\xim, \xjc}
      \frac{
         \vat{\ux}{\xic , \xjc}
         +
         \vat{\ux}{\ximm, \xjc}
      }{2}
   }{\Delta x_{\xic}} \\
   & =
   \frac{
      \vat{
         \dintrpa{\ux}{x}
      }{\xip, \xjc}
      \frac{
         \vat{\ux}{\xipp, \xjc}
         \vat{\ux}{\xic,  \xjc}
      }{2}
      -
      \vat{
         \dintrpa{\ux}{x}
      }{\xim, \xjc}
      \frac{
         \vat{\ux}{\xic,  \xjc}
         \vat{\ux}{\ximm, \xjc}
      }{2}
   }{\Delta x_{\xic}}
   +
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
   & = \dder{
      \dintrpa{\ux}{x}
      q_{xx}
   }{x}
   +
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
   }{\Delta x_{\xic}},

while the second term leads to

.. math::

   & \vat{\ux}{\xic, \xjc}
   \frac{
      \vat{
         \dintrpu{\uy}{x}
      }{\xic, \xjp}
      \frac{
         \vat{\ux}{\xic, \xjpp}
         +
         \vat{\ux}{\xic, \xjc }
      }{2}
      -
      \vat{
         \dintrpu{\uy}{x}
      }{\xic, \xjm}
      \frac{
         \vat{\ux}{\xic, \xjc }
         +
         \vat{\ux}{\xic, \xjmm}
      }{2}
   }{\Delta y} \\
   & =
   \frac{
      \vat{
         \dintrpu{\uy}{x}
      }{\xic, \xjp}
      \frac{
         \vat{\uy}{\xic, \xjpp}
         \vat{\uy}{\xic, \xjc }
      }{2}
      -
      \vat{
         \dintrpu{\uy}{x}
      }{\xic, \xjm}
      \frac{
         \vat{\uy}{\xic, \xjc }
         \vat{\uy}{\xic, \xjmm}
      }{2}
   }{\Delta y}
   +
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
   }{\Delta y} \\
   & =
   \dder{
      \dintrpu{\uy}{x}
      q_{xy}
   }{y}
   +
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
   }{\Delta y}.

Here I introduce the quadratic quantities :math:`q_{xx}` and :math:`q_{xy}`, which are defined as

.. math::

   \vat{q_{xx}}{\xip, \xjc}
   \equiv
   \frac{1}{2}
   \vat{\ux}{\xipp, \xjc}
   \vat{\ux}{\xic,  \xjc},

and

.. math::

   \vat{q_{xy}}{\xic, \xjp}
   \equiv
   \frac{1}{2}
   \vat{\ux}{\xic, \xjpp}
   \vat{\ux}{\xic, \xjc },

respectively, which are the products of the two neighbouring velocities instead of the squared velocities.

In analogy to the kinetic energy transport in the continuous domain, I request the sum of

.. math::

   \sum_{i} \sum_{j}
   \left(
      \dder{
         \dintrpa{\ux}{x}
         q_{xx}
      }{x}
      +
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
      }{\Delta x_{\xic}}
   \right)
   \vat{\Delta x}{\xic}
   \Delta y

and

.. math::

   \sum_{i} \sum_{j}
   \left(
      \dder{
         \dintrpu{\uy}{x}
         q_{xy}
      }{y}
      +
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
      }{\Delta y}
   \right)
   \vat{\Delta x}{\xic}
   \Delta y

are conserved.

The first terms are conservative and thus they inherently satisfy the requirement.
The sum of the second terms yield

.. math::

   \frac{
      \vat{\ux^2}{\xic, \xjc}
   }{2}
   \left(
      \frac{
         \vat{
            \dintrpa{\ux}{x}
         }{\xip, \xjc}
         -
         \vat{
            \dintrpa{\ux}{x}
         }{\xim, \xjc}
      }{\Delta x_{\xic}}
      +
      \frac{
         \vat{
            \dintrpu{\uy}{x}
         }{\xic, \xjp}
         -
         \vat{
            \dintrpu{\uy}{x}
         }{\xic, \xjm}
      }{\Delta y}
   \right),

which I request to vanish.

The first term in the parenthesis can be written as

.. math::

   \frac{\Delta x_{\xip}}{2 \Delta x_{\xic}} \vat{\dder{\ux}{x}}{\xip, \xjc}
   +
   \frac{\Delta x_{\xim}}{2 \Delta x_{\xic}} \vat{\dder{\ux}{x}}{\xim, \xjc}.

Although the explicit forms of :math:`\vat{\dintrpu{\uy}{x}}{\xic, \xjp}` and :math:`\vat{\dintrpu{\uy}{x}}{\xic, \xjm}` are unknown, they should be written as the linear combinations of the neighbouring velocities:

.. math::

   \vat{\dintrpu{\uy}{x}}{\xic, \cdots}
   =
   \vat{C}{\xip} \vat{\uy}{\xip, \cdots}
   +
   \vat{C}{\xim} \vat{\uy}{\xim, \cdots},

giving

.. math::

   \vat{C}{\xip} \vat{\dder{\uy}{y}}{\xip, \xjc}
   +
   \vat{C}{\xim} \vat{\dder{\uy}{y}}{\xim, \xjc}.

Thus, I notice that, they vanish if

.. math::

   \vat{C}{\xip}
   & =
   \frac{\Delta x_{\xip}}{2 \Delta x_{\xic}}, \\
   \vat{C}{\xim}
   & =
   \frac{\Delta x_{\xim}}{2 \Delta x_{\xic}}

hold because of the incompressibility constraint, and find that the placeholder should be the volume average:

.. math::

   \vat{\dintrpv{\uy}{x}}{\xic, \cdots}
   =
   \vat{C}{\xip} \vat{\uy}{\xip, \cdots}
   +
   \vat{C}{\xim} \vat{\uy}{\xim, \cdots}.

Similarly, in the :math:`y` direction at :math:`\left( \yic, \yjc \right)` where :math:`\uy` is defined, I have

.. math::

   \uy \left(
      \dder{
         \dintrpa{\ux}{y}
         \dintrpu{\uy}{x}
      }{x}
      +
      \dder{
         \dintrpa{\uy}{y}
         \dintrpa{\uy}{y}
      }{y}
   \right).

The second term is simply

.. math::

   \dder{\dintrpa{\uy}{y} q_{yy}}{y}
   +
   \frac{
      \vat{\uy^2}{\yic, \yjc}
   }{2}
   \frac{1}{2} \left(
      \vat{\dder{\uy}{y}}{\yic, \yjp}
      -
      \vat{\dder{\uy}{y}}{\yic, \yjm}
   \right),

where :math:`q_{yy}` is the quadratic quantity defined as

.. math::

   \vat{q_{yy}}{\yic, \yjp}
   \equiv
   \frac{1}{2}
   \vat{\uy}{\yic, \yjpp}
   \vat{\uy}{\yic, \yjc }.

Regarding the first term, I let the coefficients of the unknown interpolations as

.. math::

   \vat{
      \dintrpu{\uy}{x}
   }{\yim, \yjc}
   & \equiv
   \vat{c^-}{\yimm} \vat{\uy}{\yimm}
   +
   \vat{c^-}{\yic } \vat{\uy}{\yic }, \\
   \vat{
      \dintrpu{\uy}{x}
   }{\yip, \yjc}
   & \equiv
   \vat{c^+}{\yic } \vat{\uy}{\yic }
   +
   \vat{c^+}{\yipp} \vat{\uy}{\yipp},

giving

.. math::

   \color{blue}{\vat{\uy}{\yic, \yjc}}
   \frac{
      \vat{
         \dintrpa{\ux}{y}
      }{\yip, \yjc}
      \left(
         \color{blue}{
         \vat{c^+}{\yipp} \vat{\uy}{\yipp, \yjc}
         }
         +
         \vat{c^+}{\yic } \vat{\uy}{\yic , \yjc}
      \right)
      -
      \vat{
         \dintrpa{\ux}{y}
      }{\yim, \yjc}
      \left(
         \vat{c^-}{\yic } \vat{\uy}{\yic , \yjc}
         +
         \color{blue}{
         \vat{c^-}{\yimm} \vat{\uy}{\yimm, \yjc}
         }
      \right)
   }{\Delta x_{\yic}}.

I notice two constraints to identify the coefficients.

==================
Quadratic quantity
==================

I focus on the terms coloured in blue to define the quadratic quantity, which request the coefficients :math:`\vat{c^+}{\yipp}` and :math:`\vat{c^-}{\yimm}` to be :math:`1/2`, so that

.. math::

   \vat{q_{yx}}{\yip, \yjc}
   \equiv
   \frac{1}{2}
   \vat{\uy}{\yipp, \yjc}
   \vat{\uy}{\yic,  \yjc}

can be defined and I am able to make the bluish terms conservative

.. math::

   \dder{
      \dintrpa{\ux}{y}
      q_{yx}
   }{x}.

========
Residual
========

The other term yields

.. math::

   \vat{\uy^2}{\yic, \yjc}
   \frac{
      \vat{c^+}{\yic}
      \vat{
         \dintrpa{\ux}{y}
      }{\yip, \yjc}
      -
      \vat{c^-}{\yic}
      \vat{
         \dintrpa{\ux}{y}
      }{\yim, \yjc}
   }{\Delta x_{\yic}}.

To make it canceled out with the other residual, I notice the coefficients must be :math:`1/2` again.

Thus, I notice that the arithmetic average :math:`\dintrpa{q}{x}` should be used for the unknown interpolations in the :math:`y` momentum advection.

Finally, I conclude that the advective terms in the divergence form are

.. math::

   \der{
      \ux
      \ux
   }{x}
   +
   \der{
      \uy
      \ux
   }{y}
   & =
   \color{red}{
      \dder{
         \dintrpa{\ux}{x}
         \dintrpa{\ux}{x}
      }{x}
      +
      \dder{
         \dintrpv{\uy}{x}
         \dintrpa{\ux}{y}
      }{y}
   }, \\
   \der{
      \ux
      \uy
   }{x}
   +
   \der{
      \uy
      \uy
   }{y}
   & =
   \color{red}{
      \dder{
         \dintrpa{\ux}{y}
         \dintrpa{\uy}{x}
      }{x}
      +
      \dder{
         \dintrpa{\uy}{y}
         \dintrpa{\uy}{y}
      }{y}
   }.

