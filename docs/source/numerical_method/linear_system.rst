
.. _linear_system:

#############
Linear System
#############

To :ref:`treat the diffusive terms implicitly <implicit_treatment>`, we need to solve a linear system:

.. math::

    \left\{
        1
        -
        C
        \frac{1}{J}
        \dif{}{\gcs{i}}
        \left(
            \frac{J}{\sfact{i}}
            \frac{1}{\sfact{i}}
            \dif{}{\gcs{i}}
        \right)
    \right\}
    p
    =
    q,

where :math:`C` is a coefficient (e.g., coming from the diffusivity and time-step sizes).

*****************
Laplace Operators
*****************

Note that linear systems involve discrete Laplace operators:

.. math::

    \frac{1}{J}
    \dif{}{\gcs{i}}
    \left(
        \frac{J}{\sfact{i}}
        \frac{1}{\sfact{i}}
        \dif{}{\gcs{i}}
    \right),

which are pre-computed as follows.

Wall-normal velocity:

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: vector laplacian in x

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: vector laplacian in y

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: vector laplacian in z

Stream-wise velocity:

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: vector laplacian in x

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: vector laplacian in y

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: vector laplacian in z

Span-wise velocity:

.. myliteralinclude:: /../../src/fluid/predict/uz.c
    :language: c
    :tag: vector laplacian in x

.. myliteralinclude:: /../../src/fluid/predict/uz.c
    :language: c
    :tag: vector laplacian in y

.. myliteralinclude:: /../../src/fluid/predict/uz.c
    :language: c
    :tag: vector laplacian in z

Temperature:

.. myliteralinclude:: /../../src/fluid/predict/t.c
    :language: c
    :tag: scalar laplacian in x

.. myliteralinclude:: /../../src/fluid/predict/t.c
    :language: c
    :tag: scalar laplacian in y

.. myliteralinclude:: /../../src/fluid/predict/t.c
    :language: c
    :tag: scalar laplacian in z

**********************
Solving Linear Systems
**********************

Wall-normal velocity:

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: solve linear systems in x

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: solve linear systems in y

.. myliteralinclude:: /../../src/fluid/predict/ux.c
    :language: c
    :tag: solve linear systems in z

Stream-wise velocity:

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: solve linear systems in x

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: solve linear systems in y

.. myliteralinclude:: /../../src/fluid/predict/uy.c
    :language: c
    :tag: solve linear systems in z

Span-wise velocity:

.. myliteralinclude:: /../../src/fluid/predict/uz.c
    :language: c
    :tag: solve linear systems in x

.. myliteralinclude:: /../../src/fluid/predict/uz.c
    :language: c
    :tag: solve linear systems in y

.. myliteralinclude:: /../../src/fluid/predict/uz.c
    :language: c
    :tag: solve linear systems in z

Temperature:

.. myliteralinclude:: /../../src/fluid/predict/t.c
    :language: c
    :tag: solve linear systems in x

.. myliteralinclude:: /../../src/fluid/predict/t.c
    :language: c
    :tag: solve linear systems in y

.. myliteralinclude:: /../../src/fluid/predict/t.c
    :language: c
    :tag: solve linear systems in z

We need to solve :ref:`tri-diagonal matrix <tdm>`.

