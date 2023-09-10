I consider to average the above equation in the :math:`x` direction, from the left to the right walls, yielding

.. math::

   Nu &= \frac{1}{l_x l_y} \sum_{\pic} \sum_{\pjc} \left(
      \Delta x_{\pip} \Delta y_{\pjc} \sqrt{Pr} \sqrt{Ra} \, \vat{\ux}{\pip,\pjc} \vat{\dintrpa{T}{\gx}}{\pip,\pjc}
   \right) \\
   & - \frac{1}{l_x l_y} \sum_{\pjc} \left(
        \Delta y_{\pjc} \vat{T}{\text{right wall},\pjc}
      - \Delta y_{\pjc} \vat{T}{\text{left  wall},\pjc}
   \right).

Because of the normalisations, :math:`l_x` and :math:`\vat{T}{\text{left  wall},\pjc} - \vat{T}{\text{right wall},\pjc}` are unity.
Thus I have

.. math::

   Nu = \frac{1}{l_x l_y} \sum_{\pic} \sum_{\pjc} \left(
      \Delta x_{\pip} \Delta y_{\pjc} \sqrt{Pr} \sqrt{Ra} \, \vat{\ux}{\pip,\pjc} \vat{\dintrpa{T}{\gx}}{\pip,\pjc}
   \right)
   + 1,

or equivalently

.. math::

   Nu = \frac{1}{l_x l_y} \sum_{\forall \ux \text{positions}, \left( \pip, \pjc \right)} \left(
      \sqrt{Pr} \sqrt{Ra} \, J_{\pip,\pjc} \vat{\ux}{\pip,\pjc} \vat{\dintrpa{T}{\gx}}{\pip,\pjc}
   \right)
   + 1.

