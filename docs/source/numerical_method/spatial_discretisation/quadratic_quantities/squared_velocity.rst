##################################
Global Balance of Squared Velocity
##################################

******************
Sink (Dissipation)
******************

There are 9 terms dissipating the squared velocity, which are computed separately.

.. math::

    \sumzc
    \sumyc
    \sumxc
    \dkdis{1}{1}

.. myliteralinclude:: /../../src/logging/dissipated_squared_velocity.c
    :language: c
    :tag: duxdx component

.. math::

    \sumzc
    \sumyf
    \sumxf
    \dkdis{2}{1}

.. myliteralinclude:: /../../src/logging/dissipated_squared_velocity.c
    :language: c
    :tag: duxdy component

.. math::

    \sumzf
    \sumyc
    \sumxf
    \dkdis{3}{1}

.. myliteralinclude:: /../../src/logging/dissipated_squared_velocity.c
    :language: c
    :tag: duxdz component

.. math::

    \sumzc
    \sumyf
    \sumxf
    \dkdis{1}{2}

.. myliteralinclude:: /../../src/logging/dissipated_squared_velocity.c
    :language: c
    :tag: duydx component

.. math::

    \sumzc
    \sumyc
    \sumxc
    \dkdis{2}{2}

.. myliteralinclude:: /../../src/logging/dissipated_squared_velocity.c
    :language: c
    :tag: duydy component

.. math::

    \sumzf
    \sumyf
    \sumxc
    \dkdis{3}{2}

.. myliteralinclude:: /../../src/logging/dissipated_squared_velocity.c
    :language: c
    :tag: duydz component

.. math::

    \sumzf
    \sumyc
    \sumxf
    \dkdis{1}{3}

.. myliteralinclude:: /../../src/logging/dissipated_squared_velocity.c
    :language: c
    :tag: duzdx component

.. math::

    \sumzf
    \sumyf
    \sumxc
    \dkdis{2}{3}

.. myliteralinclude:: /../../src/logging/dissipated_squared_velocity.c
    :language: c
    :tag: duzdy component

.. math::

    \sumzc
    \sumyc
    \sumxc
    \dkdis{3}{3}

.. myliteralinclude:: /../../src/logging/dissipated_squared_velocity.c
    :language: c
    :tag: duzdz component

******************
Source (Injection)
******************

.. math::

    \sumzc
    \sumyc
    \sumxc
    J
    \vel{1}
    \ave{T}{\gcs{1}}

increases the total amount of the squared velocity.
Physically this term accounts for the kinetic energy injection due to the buoyancy effects.

This quantity is implemented to monitor as follows:

.. myliteralinclude:: /../../src/logging/injected_squared_velocity.c
    :language: c
    :tag: compute injected squared velocity

