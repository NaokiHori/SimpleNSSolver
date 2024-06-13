######################################
Relations involving span-wise velocity
######################################

****************
Differentiations
****************

.. math::

    \sumzf
    \sumyc
    \sumxc
    \vel{3}
    \dif{q}{\gcs{1}}
    =
    -
    \sumzf
    \sumyc
    \sumxf
    \dif{\vel{3}}{\gcs{1}}
    q,

where :math:`\vat{\vel{3}}{\frac{1}{2},\ccindex{j},\cpindex{k}} = \vat{\vel{3}}{\ngp{1} + \frac{1}{2},\ccindex{j},\cpindex{k}} = 0` is assumed (i.e., the walls do not move in the :math:`z` direction).

.. math::

    \sumzf
    \sumyc
    \sumxc
    \vel{3}
    \dif{q}{\gcs{2}}
    =
    -
    \sumzf
    \sumyf
    \sumxc
    \dif{\vel{3}}{\gcs{2}}
    q.

.. math::

    \sumzf
    \sumyc
    \sumxc
    \vel{3}
    \dif{q}{\gcs{3}}
    =
    -
    \sumzc
    \sumyc
    \sumxc
    \dif{\vel{3}}{\gcs{3}}
    q.

********
Averages
********

.. math::

    \sumzf
    \sumyc
    \sumxc
    \vel{3}
    \ave{q}{\gcs{1}}
    =
    \sumzf
    \sumyc
    \left(
        \vat{\vel{3}}{1}
        \frac{\vat{q}{\frac{1}{2}}}{2}
        +
        \sum_{i = \frac{3}{2}}^{\ngp{1} - \frac{1}{2}}
        \ave{\vel{3}}{\gcs{1}}
        q
        +
        \vat{\vel{3}}{\ngp{1}}
        \frac{\vat{q}{\ngp{1} + \frac{1}{2}}}{2}
    \right).


.. math::

    \sumzf
    \sumyc
    \sumxc
    \vel{3}
    \ave{q}{\gcs{2}}
    =
    \sumzf
    \sumyf
    \sumxc
    \ave{\vel{3}}{\gcs{2}}
    q.

.. math::

    \sumzf
    \sumyc
    \sumxc
    \vel{3}
    \ave{q}{\gcs{3}}
    =
    \sumzc
    \sumyc
    \sumxc
    \ave{\vel{3}}{\gcs{3}}
    q.

