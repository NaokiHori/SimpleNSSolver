
.. _discrete_nusselt:

##############
Nusselt number
##############

*************
Heat transfer
*************

To start, we define the heat transfer: surface-integrating :ref:`the internal energy balance <discrete_internal_energy_balance>` in the homogeneous directions yields

.. math::

    \sumzc
    \sumyc
    \frac{J}{\sfact{1}}
    \pder{T}{t}
    =
    \sumzc
    \sumyc
    \frac{1}{\sfact{1}}
    \dif{}{\gcs{1}}
    \left(
        -
        \frac{J}{\sfact{1}}
        \vel{1}
        \ave{T}{\gcs{1}}
        +
        \frac{1}{\sqrt{Pr} \sqrt{Ra}}
        \frac{J}{\sfact{1}}
        \frac{1}{\sfact{1}}
        \dif{T}{\gcs{1}}
    \right).

Assuming the flow field achieves a statistically-steady state

.. math::

    \pder{T}{t}
    \rightarrow
    0

leads to

.. math::

    \frac{1}{\sfact{1}}
    \dif{}{\gcs{1}}
    \left\{
        \sumzc
        \sumyc
        \left(
            \frac{J}{\sfact{1}}
            \vel{1}
            \ave{T}{\gcs{1}}
            -
            \frac{1}{\sqrt{Pr} \sqrt{Ra}}
            \frac{J}{\sfact{1}}
            \frac{1}{\sfact{1}}
            \dif{T}{\gcs{1}}
        \right)
    \right\}
    =
    0,

where the wall-normal differentiation and the homogeneous summations are interchanged, which is justified for rectilinear coordinates.
We introduce

.. _eq_heat_transfer:

.. math::

    \heattransfer
    \equiv
    \sumzc
    \sumyc
    \left(
        \frac{J}{\sfact{1}}
        \vel{1}
        \ave{T}{\gcs{1}}
        -
        \frac{1}{\sqrt{Pr} \sqrt{Ra}}
        \frac{J}{\sfact{1}}
        \frac{1}{\sfact{1}}
        \dif{T}{\gcs{1}}
    \right),

which is the internal energy going through a specific wall-normal position per unit time (heat transfer).
Note that the sign is decided such that *normal* cases :math:`\vat{T}{\frac{1}{2}} > \vat{T}{\ngp{1} + \frac{1}{2}}` give positive value.

The computation of the heat transfer on the walls :math:`\heattransfer` are implemented as follows:

.. myliteralinclude:: /../../src/logging/heat_transfer.c
    :language: c
    :tag: compute heat transfer on the walls

By further integrating the differential equation in the wall-normal direction, we obtain

.. math::

    \sum_{i = 1}^{\chi}
    \sfact{1}
    \frac{1}{\sfact{1}}
    \dif{\heattransfer}{\gcs{1}}
    =
    -
    \vat{\heattransfer}{\frac{1}{2}}
    +
    \vat{\heattransfer}{\xi + \frac{1}{2}}
    =
    0,

indicating that :math:`\heattransfer` is constant for all wall-normal cell faces (:math:`\frac{1}{2}, \frac{3}{2}, \cdots, \ngp{1} - \frac{1}{2}, \ngp{1} + \frac{1}{2}`).

**************
Nusselt number
**************

We focus on how much the heat transfer is enhanced compared to the reference case :math:`\heattransfer_{ref}` if the flow fields were stationary with the given :math:`Ra` and :math:`Pr`.
For stationary flow fields, heat is purely conducted diffusively (i.e., :math:`\vel{1} \equiv 0`), and the temperature profile is linear in the wall-normal direction:

.. math::

    \frac{1}{\sfact{1}}
    \dif{T}{\gcs{1}}
    =
    -
    1
    \,\,
    (\because \text{boundary conditions}),

and thus the reference heat transfer is given by

.. math::

    \heattransfer_{ref}
    =
    \sumzc
    \sumyc
    \frac{1}{\sqrt{Pr} \sqrt{Ra}}
    \frac{J}{\sfact{1}}.

The Nusselt number is defined as the ratio of them:

.. math::

    Nu
    \equiv
    \frac{\heattransfer}{\heattransfer_{ref}}.

As proved in :ref:`the global balance of squared temperature <global_balance_squared_temperature>`, :math:`\heattransfer` is linked to the source and sink of the squared temperature relation.

*******************
Squared Temperature
*******************

Now we revisit the relations derived in :ref:`the global balance of squared temperature <global_balance_squared_temperature>`:

.. math::

   \dhinjall.

Due to the Dirichlet boundary condition with respect to the temperature, we can extract :math:`T` out of summation symbols to obtain

.. math::

    \vat{
        T
    }{\frac{1}{2}}
    \vat{
        \heattransfer
    }{\frac{1}{2}}
    -
    \vat{
        T
    }{\ngp{1} + \frac{1}{2}}
    \vat{
        \heattransfer
    }{\ngp{1} + \frac{1}{2}}.

Additionally, since we fix the temperature difference of the two walls

.. math::

    \vat{T}{\frac{1}{2}}
    -
    \vat{T}{\ngp{1} + \frac{1}{2}}
    \equiv
    1,

and :math:`\heattransfer` is equal at every wall-normal cell faces

.. math::

    \heattransfer
    \equiv
    \vat{\heattransfer}{\frac{1}{2}}
    =
    \vat{\heattransfer}{\ngp{1} + \frac{1}{2}},

we notice

.. math::

    \heattransfer
    =
    \dhinjall,

and of course

.. math::

    \heattransfer
    =
    \dhdisall.

****************
Squared Velocity
****************

Integrating :ref:`the definition of heat transfer <eq_heat_transfer>` in the wall-normal direction yields

.. math::

    \sumxc
    \sfact{1}
    \heattransfer
    =
    \sumzc
    \sumyc
    \sumxc
    J
    \vel{1}
    \ave{T}{\gcs{1}}
    -
    \sumzc
    \sumyc
    \sumxc
    \frac{1}{\sqrt{Pr} \sqrt{Ra}}
    \frac{J}{\sfact{1}}
    \dif{T}{\gcs{1}}.

The left-hand side is equal to :math:`\heattransfer` since :math:`\sumxc \sfact{1} \equiv 1` (recall that we assume the wall-normal length of the domain is unity).
The second term in the right-hand side leads to

.. math::

    &
    -
    \sumzc
    \sumyc
    \frac{1}{\sqrt{Pr} \sqrt{Ra}}
    \frac{J}{\sfact{1}}
    \sumxc
    \dif{T}{\gcs{1}}

    =
    &
    -
    \sumzc
    \sumyc
    \frac{1}{\sqrt{Pr} \sqrt{Ra}}
    \frac{J}{\sfact{1}}
    \left(
        \vat{T}{\ngp{1} + \frac{1}{2}}
        -
        \vat{T}{\frac{1}{2}}
    \right)

    =
    &
    \sumzc
    \sumyc
    \frac{1}{\sqrt{Pr} \sqrt{Ra}}
    \frac{J}{\sfact{1}}

    =
    &
    \heattransfer_{ref},

since :math:`J / \sfact{1}` is independent to the homogeneous directions.

Thus we notice that

.. math::

    \heattransfer
    =
    \sumzc
    \sumyc
    \sumxc
    J
    \vel{1}
    \ave{T}{\gcs{1}}
    +
    \heattransfer_{ref},

which relates the Nusselt number with the squared velocity relations.

