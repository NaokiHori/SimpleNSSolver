
.. _quadratic_quantities:

####################
Quadratic Quantities
####################

Here we aim at confirming if :ref:`the theoretical relations of the quadratic quantities <quadratic_quantity_balance>` are satisfied even after discretized.

Specifically, as discrete counterparts, we consider :ref:`the discrete momentum balance <discrete_momentum_balance>` multiplied by the Jacobian determinant (incorporating the volume integrals) and the velocity:

.. math::

    \sumzc
    \sumyc
    \sumxf
    J
    \vel{1}
    \pder{\vel{1}}{t}
    =
    &
    -
    \sumzc
    \sumyc
    \sumxf
    J
    \vel{1}
    \dmomadv{1}{1}
    -
    \sumzc
    \sumyc
    \sumxf
    J
    \vel{1}
    \dmomadv{2}{1}
    -
    \sumzc
    \sumyc
    \sumxf
    J
    \vel{1}
    \dmomadv{3}{1}

    &
    -
    \sumzc
    \sumyc
    \sumxf
    J
    \vel{1}
    \dmompre{1}

    &
    +
    \sumzc
    \sumyc
    \sumxf
    J
    \vel{1}
    \dmomdif{1}{1}
    +
    \sumzc
    \sumyc
    \sumxf
    J
    \vel{1}
    \dmomdif{2}{1}
    +
    \sumzc
    \sumyc
    \sumxf
    J
    \vel{1}
    \dmomdif{3}{1}

    &
    +
    \sumzc
    \sumyc
    \sumxf
    J
    \vel{1}
    \ave{T}{\gcs{1}},

.. math::

    \sumzc
    \sumyf
    \sumxc
    J
    \vel{2}
    \pder{\vel{2}}{t}
    =
    &
    -
    \sumzc
    \sumyf
    \sumxc
    J
    \vel{2}
    \dmomadv{1}{2}
    -
    \sumzc
    \sumyf
    \sumxc
    J
    \vel{2}
    \dmomadv{2}{2}
    -
    \sumzc
    \sumyf
    \sumxc
    J
    \vel{2}
    \dmomadv{3}{2}

    &
    -
    \sumzc
    \sumyf
    \sumxc
    J
    \vel{2}
    \dmompre{2}

    &
    +
    \sumzc
    \sumyf
    \sumxc
    J
    \vel{2}
    \dmomdif{1}{2}
    +
    \sumzc
    \sumyf
    \sumxc
    J
    \vel{2}
    \dmomdif{2}{2}
    +
    \sumzc
    \sumyf
    \sumxc
    J
    \vel{2}
    \dmomdif{3}{2},

.. math::

    \sumzf
    \sumyc
    \sumxc
    J
    \vel{3}
    \pder{\vel{3}}{t}
    =
    &
    -
    \sumzf
    \sumyc
    \sumxc
    J
    \vel{3}
    \dmomadv{1}{3}
    -
    \sumzf
    \sumyc
    \sumxc
    J
    \vel{3}
    \dmomadv{2}{3}
    -
    \sumzf
    \sumyc
    \sumxc
    J
    \vel{3}
    \dmomadv{3}{3}

    &
    -
    \sumzf
    \sumyc
    \sumxc
    J
    \vel{3}
    \dmompre{3}

    &
    +
    \sumzf
    \sumyc
    \sumxc
    J
    \vel{3}
    \dmomdif{1}{3}
    +
    \sumzf
    \sumyc
    \sumxc
    J
    \vel{3}
    \dmomdif{2}{3}
    +
    \sumzf
    \sumyc
    \sumxc
    J
    \vel{3}
    \dmomdif{3}{3}.

Similarly, multiplying :ref:`the discrete internal energy balance <discrete_internal_energy_balance>` by temperature and volume-integrating it yield

.. math::

    \sumzc
    \sumyc
    \sumxc
    J
    T
    \pder{T}{t}
    =
    &
    -
    \sumzc
    \sumyc
    \sumxc
    J
    T
    \dtempadv{1}
    -
    \sumzc
    \sumyc
    \sumxc
    J
    T
    \dtempadv{2}
    -
    \sumzc
    \sumyc
    \sumxc
    J
    T
    \dtempadv{3}

    &
    +
    \sumzc
    \sumyc
    \sumxc
    J
    T
    \dtempdif{1}
    +
    \sumzc
    \sumyc
    \sumxc
    J
    T
    \dtempdif{2}
    +
    \sumzc
    \sumyc
    \sumxc
    J
    T
    \dtempdif{3}.

It should be noted that, since the velocities are defined at different positions due to :ref:`the staggered grid arrangement <domain_setup>`, quadratic quantities are also evaluated at different locations, which is the reason why we consider three components separately.

With some algebra, we obtain

.. math::

    &
    \sumzc
    \sumyc
    \sumxf
    J
    \pder{k_1}{t}
    +
    \sumzc
    \sumyf
    \sumxc
    J
    \pder{k_2}{t}
    +
    \sumzf
    \sumyc
    \sumxc
    J
    \pder{k_3}{t}

    =
    &
    -
    \sumzc
    \sumyc
    \sumxc
    \dkdis{1}{1}
    -
    \sumzc
    \sumyf
    \sumxf
    \dkdis{2}{1}
    -
    \sumzf
    \sumyc
    \sumxf
    \dkdis{3}{1}

    &
    -
    \sumzc
    \sumyf
    \sumxf
    \dkdis{1}{2}
    -
    \sumzc
    \sumyc
    \sumxc
    \dkdis{2}{2}
    -
    \sumzf
    \sumyf
    \sumxc
    \dkdis{3}{2}

    &
    -
    \sumzf
    \sumyc
    \sumxf
    \dkdis{1}{3}
    -
    \sumzf
    \sumyf
    \sumxc
    \dkdis{2}{3}
    -
    \sumzc
    \sumyc
    \sumxc
    \dkdis{3}{3}

    &
    +
    \sumzc
    \sumyc
    \sumxf
    J
    \vel{1}
    \ave{T}{\gcs{1}},

and

.. math::

    &
    \sumzc
    \sumyc
    \sumxc
    J
    \pder{h}{t}

    =
    &
    -
    \sumzc
    \sumyc
    \sumxf
    \dhdis{1}
    -
    \sumzc
    \sumyf
    \sumxc
    \dhdis{2}
    -
    \sumzf
    \sumyc
    \sumxc
    \dhdis{3}

    &
    -
    \sumzc
    \sumyc
    \frac{1}{\sqrt{Pr} \sqrt{Ra}}
    \vat{
        \left(
            \frac{J}{\sfact{1}}
            T
            \frac{1}{\sfact{1}}
            \dif{T}{\gcs{1}}
        \right)
    }{\frac{1}{2}}

    &
    +
    \sumzc
    \sumyc
    \frac{1}{\sqrt{Pr} \sqrt{Ra}}
    \vat{
        \left(
            \frac{J}{\sfact{1}}
            T
            \frac{1}{\sfact{1}}
            \dif{T}{\gcs{1}}
        \right)
    }{\ngp{1} + \frac{1}{2}},

with respect to the squared velocity and the squared temperature, respectively.

To derive them, the individual components are focused below.

.. toctree::
    :maxdepth: 1

    derivation/prerequisite/main
    derivation/advective
    derivation/pressure_gradient
    derivation/diffusive

The derived two relations have source terms and sink terms, which are elaborated below as well as their implementations.

.. toctree::
    :maxdepth: 1

    squared_velocity
    squared_temperature

