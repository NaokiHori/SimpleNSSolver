###############
Diffusive terms
###############

We utilise :ref:`the relations listed in the prerequisite <energy_prerequisite>` to derive the following relations.

****************
Squared velocity
****************

From the wall-normal momentum balance, we obtain

.. math::

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
    \dkdis{3}{1}.

From the stream-wise momentum balance, we obtain

.. math::

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
    \dkdis{3}{2}.

From the span-wise momentum balance, we obtain

.. math::

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
    \dkdis{3}{3}.

In total, all terms work to dissipate :math:`k`.

*******************
Squared temperature
*******************

We obtain

.. math::

    -
    \sumzc
    \sumyc
    \sumxf
    J
    \frac{1}{\sqrt{Pr} \sqrt{Ra}}
    \left(
        \frac{1}{\sfact{1}}
        \dif{T}{\gcs{1}}
    \right)^2
    -
    \sumzc
    \sumyf
    \sumxc
    J
    \frac{1}{\sqrt{Pr} \sqrt{Ra}}
    \left(
        \frac{1}{\sfact{2}}
        \dif{T}{\gcs{2}}
    \right)^2
    -
    \sumzf
    \sumyc
    \sumxc
    J
    \frac{1}{\sqrt{Pr} \sqrt{Ra}}
    \left(
        \frac{1}{\sfact{3}}
        \dif{T}{\gcs{3}}
    \right)^2

as the dissipative terms, while

.. math::

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
    }{\ngp{1} + \frac{1}{2}}

as the conduction on the walls.

