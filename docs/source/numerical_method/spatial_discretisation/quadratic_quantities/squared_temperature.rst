
.. _global_balance_squared_temperature:

#####################################
Global Balance of Squared Temperature
#####################################

******************
Sink (Dissipation)
******************

There are 3 terms dissipating the squared temperature, which are computed separately.

.. math::

    \sumzc
    \sumyc
    \sumxf
    \dhdis{1}

.. myliteralinclude:: /../../src/logging/dissipated_squared_temperature.c
    :language: c
    :tag: dtdx component

.. math::

    \sumzc
    \sumyf
    \sumxc
    \dhdis{2}

.. myliteralinclude:: /../../src/logging/dissipated_squared_temperature.c
    :language: c
    :tag: dtdy component

.. math::

    \sumzf
    \sumyc
    \sumxc
    \dhdis{3}

.. myliteralinclude:: /../../src/logging/dissipated_squared_temperature.c
    :language: c
    :tag: dtdz component

******************
Source (Injection)
******************

.. math::

    \dhinjall

injects squared temperature, which will be elaborated :ref:`later <discrete_nusselt>`.

