#########
Advection
#########

********************
Wall-normal relation
********************

.. math::

    -
    \dmomadv{1}{1}

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: advected in x

.. math::

    -
    \dmomadv{2}{1}

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: advected in y

.. math::

    -
    \dmomadv{3}{1}

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: advected in z

********************
Stream-wise relation
********************

.. math::

    -
    \dmomadv{1}{2}

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: advected in x

.. math::

    -
    \dmomadv{2}{2}

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: advected in y

.. math::

    -
    \dmomadv{3}{2}

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: advected in z

******************
Span-wise relation
******************

.. math::

    -
    \dmomadv{1}{3}

.. myliteralinclude:: /../../src/fluid/predict/uz.c
    :language: c
    :tag: advected in x

.. math::

    -
    \dmomadv{2}{3}

.. myliteralinclude:: /../../src/fluid/predict/uz.c
    :language: c
    :tag: advected in y

.. math::

    -
    \dmomadv{3}{3}

.. myliteralinclude:: /../../src/fluid/predict/uz.c
    :language: c
    :tag: advected in z

