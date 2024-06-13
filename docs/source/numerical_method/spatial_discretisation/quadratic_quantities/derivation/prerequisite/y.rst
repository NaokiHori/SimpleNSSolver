########################################
Relations involving stream-wise velocity
########################################

****************
Differentiations
****************

.. math::

  \sumzc
  \sumyf
  \sumxc
  \vel{2}
  \dif{q}{\gcs{1}}
  =
  -
  \sumzc
  \sumyf
  \left(
    \vat{\left(\vel{2} q\right)}{\frac{1}{2}}
    +
    \sumxf
    \dif{\vel{2}}{\gcs{1}}
    q
    -
    \vat{\left(\vel{2} q\right)}{\ngp{1} + \frac{1}{2}}
  \right).

.. math::

  \sumzc
  \sumyf
  \sumxc
  \vel{2}
  \dif{q}{\gcs{2}}
  =
  -
  \sumzc
  \sumyc
  \sumxc
  \dif{\vel{2}}{\gcs{2}}
  q.

.. math::

  \sumzc
  \sumyf
  \sumxc
  \vel{2}
  \dif{q}{\gcs{3}}
  =
  -
  \sumzf
  \sumyf
  \sumxc
  \dif{\vel{2}}{\gcs{3}}
  q.

********
Averages
********

.. math::

  \sumzc
  \sumyf
  \sumxc
  \vel{2}
  \ave{q}{\gcs{1}}
  =
  \sumzc
  \sumyf
  \left(
    \vat{\vel{2}}{1}
    \frac{\vat{q}{\frac{1}{2}}}{2}
    +
    \sum_{i = \frac{3}{2}}^{\ngp{1} - \frac{1}{2}}
    \ave{\vel{2}}{\gcs{1}}
    q
    +
    \vat{\vel{2}}{\ngp{1}}
    \frac{\vat{q}{\ngp{1} + \frac{1}{2}}}{2}
  \right).

.. math::

  \sumzc
  \sumyf
  \sumxc
  \vel{2}
  \ave{q}{\gcs{2}}
  =
  \sumzc
  \sumyc
  \sumxc
  \ave{\vel{2}}{\gcs{2}}
  q.

.. math::

  \sumzc
  \sumyf
  \sumxc
  \vel{2}
  \ave{q}{\gcs{3}}
  =
  \sumzf
  \sumyf
  \sumxc
  \ave{\vel{2}}{\gcs{3}}
  q.

