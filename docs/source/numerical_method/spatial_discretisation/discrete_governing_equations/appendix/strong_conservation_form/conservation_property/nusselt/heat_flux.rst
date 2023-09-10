I consider to average the equation of the internal energy balance:

.. math::

   \der{\vat{T}{\pic,\pjc}}{t}
   & + \frac{
       \vat{\Delta y \ux}{\pip,\pjc} \vat{\dintrpa{T}{\gx}}{\pip,\pjc}
     - \vat{\Delta y \ux}{\pim,\pjc} \vat{\dintrpa{T}{\gx}}{\pim,\pjc}
   }{\Delta x_{\pic} \Delta y_{\pjc}} \\
   & + \frac{
       \vat{\Delta x \uy}{\pic,\pjp} \vat{\dintrpa{T}{\gy}}{\pic,\pjp}
     - \vat{\Delta x \uy}{\pic,\pjm} \vat{\dintrpa{T}{\gy}}{\pic,\pjm}
   }{\Delta x_{\pic} \Delta y_{\pjc}} \\
   & - \frac{1}{\sqrt{Pr} \sqrt{Ra}} \frac{
       \vat{\frac{\Delta y}{\Delta x} \diffe{T}{\gx}}{\pip,\pjc}
     - \vat{\frac{\Delta y}{\Delta x} \diffe{T}{\gx}}{\pim,\pjc}
   }{\Delta x_{\pic} \Delta y_{\pjc}} \\
   & - \frac{1}{\sqrt{Pr} \sqrt{Ra}} \frac{
       \vat{\frac{\Delta x}{\Delta y} \diffe{T}{\gy}}{\pic,\pjp}
     - \vat{\frac{\Delta x}{\Delta y} \diffe{T}{\gy}}{\pic,\pjm}
   }{\Delta x_{\pic} \Delta y_{\pjc}} = 0

in the :math:`y` direction and in time.

In the original coordinate system :math:`\left( x, y \right)`, integrating a quantity :math:`q` in :math:`y` direction is written as

.. math::

   \int q dy,

which is

.. math::

   \int q \der{y}{\gy} d\gy

in the general coordinate system :math:`\left( \gx, \gy \right)`.
When discretised, it is

.. math::

   \sum_j q_j \Delta y_j.

Based on this relation, I consider to average the equation of the internal energy.

The first term vanishes since I am interested in the statistically-steady state.

The third and fifth terms lead to

.. math::

   \sum_j \left(
       \vat{\uy}{\pic,\pjp} \vat{\dintrpa{T}{\gy}}{\pic,\pjp}
     - \vat{\uy}{\pic,\pjm} \vat{\dintrpa{T}{\gy}}{\pic,\pjm}
   \right)

and

.. math::

   \sum_j - \frac{1}{\sqrt{Pr} \sqrt{Ra}} \left(
       \vat{\frac{1}{\Delta y} \diffe{T}{\gy}}{\pic,\pjp}
     - \vat{\frac{1}{\Delta y} \diffe{T}{\gy}}{\pic,\pjm}
   \right),

which are zero because of the periodicity in the :math:`y` direction.

Thus I have

.. math::

   & \frac{1}{l_y} \sum_{\pjc} \frac{\Delta y_{\pjc}}{\Delta x_{\pic}} \left\{
      \left(
          \vat{\ux}{\pip,\pjc} \vat{\dintrpa{T}{\gx}}{\pip,\pjc}
        - \vat{\ux}{\pim,\pjc} \vat{\dintrpa{T}{\gx}}{\pim,\pjc}
      \right)
   \right\} \\
   & = \frac{1}{l_y} \sum_{\pjc} \frac{\Delta y_{\pjc}}{\Delta x_{\pic}} \left\{
      \frac{1}{\sqrt{Pr} \sqrt{Ra}}
      \left(
          \vat{\frac{1}{\Delta x} \diffe{T}{\gx}}{\pip,\pjc}
        - \vat{\frac{1}{\Delta x} \diffe{T}{\gx}}{\pim,\pjc}
      \right)
   \right\}.

Now I integrate this relation in the :math:`x` direction (:math:`\int q dx \approx \sum_i q_i \Delta x_i`) to a specific position :math:`i`, yielding

.. math::

   & \frac{1}{l_y} \sum_{\pic = 1}^{I} \sum_{\pjc} \Delta y_{\pjc} \left\{
      \sqrt{Pr} \sqrt{Ra} \left(
          \vat{\ux}{\pip,\pjc} \vat{\dintrpa{T}{\gx}}{\pip,\pjc}
        - \vat{\ux}{\pim,\pjc} \vat{\dintrpa{T}{\gx}}{\pim,\pjc}
      \right)
      - \left(
          \vat{\frac{1}{\Delta x} \diffe{T}{\gx}}{\pip,\pjc}
        - \vat{\frac{1}{\Delta x} \diffe{T}{\gx}}{\pim,\pjc}
      \right)
   \right\} \\
   & =
   \frac{1}{l_y} \sum_{\pjc} \Delta y_{\pjc} \left\{
      \sqrt{Pr} \sqrt{Ra} \left(
          \vat{\ux}{\pip        ,\pjc} \vat{\dintrpa{T}{\gx}}{\pip        ,\pjc}
        - \vat{\ux}{\frac{1}{2},\pjc} \vat{\dintrpa{T}{\gx}}{\frac{1}{2},\pjc}
      \right)
      - \left(
          \vat{\frac{1}{\Delta x} \diffe{T}{\gx}}{\pip        ,\pjc}
        - \vat{\frac{1}{\Delta x} \diffe{T}{\gx}}{\frac{1}{2},\pjc}
      \right)
   \right\} = 0,

i.e.

.. math::

   \frac{1}{l_y} \sum_{\pjc} \Delta y_{\pjc} \left(
      \sqrt{Pr} \sqrt{Ra} \, \vat{u}{\pip,\pjc} \vat{\dintrpa{T}{\gx}}{\pip,\pjc}
      - \vat{\frac{1}{\Delta x} \diffe{T}{\gx}}{\pip,\pjc}
   \right)
   =
   \frac{1}{l_y} \sum_{\pjc} \Delta y_{\pjc} \left(
     \sqrt{Pr} \sqrt{Ra} \, \vat{u}{\frac{1}{2},\pjc} \vat{\dintrpa{T}{\gx}}{\frac{1}{2},\pjc}
      - \vat{\frac{1}{\Delta x} \diffe{T}{\gx}}{\frac{1}{2},\pjc}
   \right).

Note that this relation holds for all :math:`i`.

The left-hand side denotes the value at a specific wall-normal location (:math:`x` cell faces where :math:`u_x` are defined), while the right-hand side describes the value on the left wall.
Since I assume the walls are impermeable, the first term (advection of the temperature) vanishes and thus I finally notice

.. math::

   Nu
   & =
   \frac{1}{l_y} \sum_{\pjc} \Delta y_{\pjc} \left(
      - \vat{\frac{1}{\Delta x} \diffe{T}{\gx}}{\text{wall},\pjc}
   \right) \\
   & =
   \frac{1}{l_y} \sum_{\pjc} \Delta y_{\pjc} \left(
      \sqrt{Pr} \sqrt{Ra} \, \vat{\ux}{\pip,\pjc} \vat{\dintrpa{T}{\gx}}{\pip,\pjc}
      - \vat{\frac{1}{\Delta x} \diffe{T}{\gx}}{\pip,\pjc}
   \right).

.. note::

   In the continuous level, this Nusselt constancy in the wall-normal direction holds for all :math:`x`.
   After discretised, this relation is satisfied only where :math:`u_x` are defined.

