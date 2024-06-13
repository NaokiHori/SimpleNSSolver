#########
Diffusion
#########

********************
Wall-normal relation
********************

.. math::

    \dmomdif{1}{1}

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: vector laplacian in x

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: diffused in x

.. math::

    \dmomdif{2}{1}

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: vector laplacian in y

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: diffused in y

.. math::

    \dmomdif{3}{1}

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: vector laplacian in z

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: diffused in z

********************
Stream-wise relation
********************

.. math::

    \dmomdif{1}{2}

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: vector laplacian in x

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: diffused in x

.. math::

    \dmomdif{2}{2}

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: vector laplacian in y

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: diffused in y

.. math::

    \dmomdif{3}{2}

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: vector laplacian in z

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: diffused in z

******************
Span-wise relation
******************

.. math::

    \dmomdif{1}{3}

.. myliteralinclude:: /../../src/fluid/predict/uz.c
    :language: c
    :tag: diffused in x

.. myliteralinclude:: /../../src/fluid/predict/uz.c
    :language: c
    :tag: vector laplacian in x

.. math::

    \dmomdif{2}{3}

.. myliteralinclude:: /../../src/fluid/predict/uz.c
    :language: c
    :tag: diffused in y

.. myliteralinclude:: /../../src/fluid/predict/uz.c
    :language: c
    :tag: vector laplacian in y

.. math::

    \dmomdif{3}{3}

.. myliteralinclude:: /../../src/fluid/predict/uz.c
    :language: c
    :tag: vector laplacian in z

.. myliteralinclude:: /../../src/fluid/predict/uz.c
    :language: c
    :tag: diffused in z

