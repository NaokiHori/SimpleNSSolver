
.. _momentum_advective_terms_divergence_form:

###############
Divergence form
###############

*******
Summary
*******

The discrete advective terms are given by

.. math::

   \dder{
      \dintrpa{\ux}{x}
      \dintrpa{\ux}{x}
   }{x}
   +
   \dder{
      \dintrpu{\uy}{x}
      \dintrpa{\ux}{y}
   }{y}

at :math:`\left( \xic, \xjc \right)`, and

.. math::

   \dder{
      \dintrpa{\ux}{y}
      \dintrpu{\uy}{x}
   }{x}
   +
   \dder{
      \dintrpa{\uy}{y}
      \dintrpa{\uy}{y}
   }{y}

at :math:`\left( \yic, \yjc \right)`, respectively.
Note that the interpolations :math:`\dintrpu{\uy}{x}` and :math:`\dintrpu{\ux}{y}` are undefined here.

**********
Derivation
**********

To begin with, I focus on the advective terms in the momentum balance

.. math::

   u_j \der{u_i}{x_j},

which is known as the gradient form and is not conservative itself.
To make this term conservative, I use the relation:

.. math::

   u_j \der{u_i}{x_j}
   +
   u_i \der{u_j}{x_j}
   =
   \der{u_j u_i}{x_j}.

Note that the second term on the left-hand side is the incompressibility constraint weighted by :math:`u_i`.
Thus, when the incompressibility constraint is satisfied, I can write the advective term in a conservative form, which is known as the divergence form.

Since they are inherently momentum-conservative, I first discretise this term:

.. math::

   \dder{
      \dintrpu{\ux}{x}
      \dintrpu{\ux}{x}
   }{x}
   +
   \dder{
      \dintrpu{\uy}{x}
      \dintrpu{\ux}{y}
   }{y}

at :math:`\left( \xic, \xjc \right)` and

.. math::

   \dder{
      \dintrpu{\ux}{y}
      \dintrpu{\uy}{x}
   }{x}
   +
   \dder{
      \dintrpu{\uy}{y}
      \dintrpu{\uy}{y}
   }{y}

at :math:`\left( \yic, \yjc \right)`, which are clearly conservative.

.. note::

   There are two possibilities to write :math:`\ux \ux`.

      * :math:`\dintrpu{\left( \ux \right)^2}{x}`: interpolated after squared.

      * :math:`\left( \dintrpu{\ux}{x} \right)^2 = \dintrpu{\ux}{x} \dintrpu{\ux}{x}`: squared after interpolated.

   I take the second option since it is consistent with the other direction :math:`\dintrpu{\uy}{x} \dintrpu{\ux}{y}`.

Placeholders :math:`\dintrpu{q}{x}` and :math:`\dintrpu{q}{y}` can be partially replaced by the arithmetic averages :math:`\dintrpa{q}{x}` and :math:`\dintrpa{q}{y}`, because cell centers are positioned in the middle of the surrounding cell faces in this project (see :ref:`the domain setup <domain_setup>`), giving

.. math::

   \dder{
      \dintrpa{\ux}{x}
      \dintrpa{\ux}{x}
   }{x}
   +
   \dder{
      \dintrpu{\uy}{x}
      \dintrpa{\ux}{y}
   }{y}

at :math:`\left( \xic, \xjc \right)` and

.. math::

   \dder{
      \dintrpa{\ux}{y}
      \dintrpu{\uy}{x}
   }{x}
   +
   \dder{
      \dintrpa{\uy}{y}
      \dintrpa{\uy}{y}
   }{y}

at :math:`\left( \yic, \yjc \right)`.

