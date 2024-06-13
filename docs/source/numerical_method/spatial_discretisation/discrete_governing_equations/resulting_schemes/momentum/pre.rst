#################
Pressure-gradient
#################

********************
Wall-normal relation
********************

.. math::

    -
    \dmompre{1}

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: pressure gradient effect

********************
Stream-wise relation
********************

.. math::

    -
    \dmompre{2}

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: pressure gradient effect

******************
Span-wise relation
******************

.. math::

    -
    \dmompre{3}

.. myliteralinclude:: /../../src/fluid/predict/uz.c
    :language: c
    :tag: pressure gradient effect

