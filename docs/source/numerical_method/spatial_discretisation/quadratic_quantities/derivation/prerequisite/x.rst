########################################
Relations involving wall-normal velocity
########################################

****************
Differentiations
****************

.. math::

    \sumzc
    \sumyc
    \sumxf
    \vel{1}
    \dif{q}{\gcs{1}}
    =
    -
    \sumzc
    \sumyc
    \sumxc
    \dif{\vel{1}}{\gcs{1}}
    q,

where :math:`\vat{\vel{1}}{\frac{1}{2}} = \vat{\vel{1}}{\ngp{1} + \frac{1}{2}} = 0` is used.

.. math::

    \sumzc
    \sumyc
    \sumxf
    \vel{1}
    \dif{q}{\gcs{2}}
    =
    -
    \sumzc
    \sumyf
    \sumxf
    \dif{\vel{1}}{\gcs{2}}
    q.

.. math::

    \sumzc
    \sumyc
    \sumxf
    \vel{1}
    \dif{q}{\gcs{3}}
    =
    -
    \sumzf
    \sumyc
    \sumxf
    \dif{\vel{1}}{\gcs{3}}
    q.

********
Averages
********

.. math::

    \sumzc
    \sumyc
    \sumxf
    \vel{1}
    \ave{q}{\gcs{1}}
    =
    \sumzc
    \sumyc
    \sumxc
    \ave{\vel{1}}{\gcs{1}}
    q.

.. math::

    \sumzc
    \sumyc
    \sumxf
    \vel{1}
    \ave{q}{\gcs{2}}
    =
    \sumzc
    \sumyf
    \sumxf
    \ave{\vel{1}}{\gcs{2}}
    q.

.. math::

    \sumzc
    \sumyc
    \sumxf
    \vel{1}
    \ave{q}{\gcs{3}}
    =
    \sumzf
    \sumyc
    \sumxf
    \ave{\vel{1}}{\gcs{3}}
    q.

