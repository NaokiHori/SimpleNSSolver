
.. include:: /references.txt

#################
Resulting schemes
#################

We discretize the mentioned equations as follows.
Note that we define several symbols for notational convenience.

.. toctree::
    :maxdepth: 1

    symbols/average
    symbols/differentiation
    symbols/summation

The advective terms are implemented in tri-diagonal-like ways to insist the advective operators are skew-symmetric (c.f., |VERSTAPPEN2003|).
Since the diffusive terms are given by the vector Laplacian of the velocity vector, the discrete Laplace operators are pre-computed.

.. _discrete_incompressibility:

****************************
Incompressibility constraint
****************************

We define this relation for each cell center :math:`\left( \ccidx{i}, \ccidx{j}, \ccidx{k} \right)`:

.. math::

    \ddiv{1}
    +
    \ddiv{2}
    +
    \ddiv{3}
    =
    0.

The left-hand side is monitored to confirm if the maximum divergence of the flow field is sufficiently small (:math:`\approx 0`):

.. myliteralinclude:: /../../src/logging/max_divergence.c
    :language: c
    :tag: check max local divergence

.. _discrete_momentum_balance:

****************
Momentum balance
****************

The following relations are defined at the corresponding cell faces.

.. math::

    \pder{\vel{1}}{t}
    =
    &
    -
    \dmomadv{1}{1}
    -
    \dmomadv{2}{1}
    -
    \dmomadv{3}{1}

    &
    -
    \dmompre{1}

    &
    +
    \dmomdif{1}{1}
    +
    \dmomdif{2}{1}
    +
    \dmomdif{3}{1}

    &
    +
    \dmombuo.

.. math::

    \pder{\vel{2}}{t}
    =
    &
    -
    \dmomadv{1}{2}
    -
    \dmomadv{2}{2}
    -
    \dmomadv{3}{2}

    &
    -
    \dmompre{2}

    &
    +
    \dmomdif{1}{2}
    +
    \dmomdif{2}{2}
    +
    \dmomdif{3}{2}.

.. math::

    \pder{\vel{3}}{t}
    =
    &
    -
    \dmomadv{1}{3}
    -
    \dmomadv{2}{3}
    -
    \dmomadv{3}{3}

    &
    -
    \dmompre{3}

    &
    +
    \dmomdif{1}{3}
    +
    \dmomdif{2}{3}
    +
    \dmomdif{3}{3}.

The implementations are elaborated below.

.. toctree::
    :maxdepth: 1

    momentum/adv
    momentum/pre
    momentum/dif
    momentum/buo

.. _discrete_internal_energy_balance:

***********************
Internal energy balance
***********************

We define this relation for each cell center :math:`\left( \ccidx{i}, \ccidx{j}, \ccidx{k} \right)`:

.. math::

    \pder{T}{t}
    =
    &
    -
    \dtempadv{1}
    -
    \dtempadv{2}
    -
    \dtempadv{3}

    &
    +
    \dtempdif{1}
    +
    \dtempdif{2}
    +
    \dtempdif{3}.

The implementations are elaborated below.

.. toctree::
    :maxdepth: 1

    t/adv
    t/dif

