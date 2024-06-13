###############
Advective terms
###############

****************
Squared velocity
****************

To begin with, we consider

.. math::

    \newcommand{\tmp}[1]{
        \dmomadv{#1}{1}
        \dmomadv{#1}{2}
        \dmomadv{#1}{3}
        -
        \frac{1}{J}
        \vel{#1}
        \ave{
            \dif{
                \left(
                    \frac{J}{\sfact{1}}
                    \vel{1}
                \right)
            }{\gcs{1}}
            +
            \dif{
                \left(
                    \frac{J}{\sfact{2}}
                    \vel{2}
                \right)
            }{\gcs{2}}
            +
            \dif{
                \left(
                    \frac{J}{\sfact{3}}
                    \vel{3}
                \right)
            }{\gcs{3}}
        }{\gcs{#1}}
    }
    \tmp{1},

    \tmp{2},

    \tmp{3},

where the first three terms for each direction are the advective terms in the momentum balance, while the last terms are the product of the corresponding velocity and the averaged incompressibility constraint and thus zero.
With some algebra, we find that they are equal to

.. math::

    \newcommand{\tmp}[1]{
        -
        \frac{1}{J}
        \dif{
           \left(
              \ave{
                 \frac{J}{\sfact{#1}} \vel{1}
              }{\gcs{#1}}
              \ave{\vel{#1}}{\gcs{1}}
           \right)
        }{\gcs{1}}
        -
        \frac{1}{J}
        \dif{
           \left(
              \ave{
                 \frac{J}{\sfact{#1}} \vel{2}
              }{\gcs{#1}}
              \ave{\vel{#1}}{\gcs{2}}
           \right)
        }{\gcs{2}}
        -
        \frac{1}{J}
        \dif{
           \left(
              \ave{
                 \frac{J}{\sfact{#1}} \vel{3}
              }{\gcs{#1}}
              \ave{\vel{#1}}{\gcs{3}}
           \right)
        }{\gcs{3}}
    }
    \tmp{1},

    \tmp{2},

    \tmp{3},

respectively.

In short, the above relation is the discrete description of

.. math::

    -
    \vel{j} \pder{\vel{i}}{x_j}
    -
    \vel{i} \pder{\vel{j}}{x_j}
    =
    -
    \pder{\vel{j} \vel{i}}{x_j},

i.e., the advective and the divergence forms of the advective terms are equivalent.

Now the global energy balance attributed to the advective terms are considered:

.. math::

    \newcommand{\tmp}[1]{
        J
        \vel{#1}
        \left(
            \dmomadv{#1}{1}
            \dmomadv{#1}{2}
            \dmomadv{#1}{3}
        \right)
    }
    \sumzc
    \sumyc
    \sumxf
    \tmp{1},

    \sumzc
    \sumyf
    \sumxc
    \tmp{2},

    \sumzf
    \sumyc
    \sumxc
    \tmp{3}.

Our objective is to prove they are all zero, indicating that the advective terms do not alter the net amount of the quadratic quantities.

By using the relations derived in :ref:`the prerequisite <energy_prerequisite>`, each term leads to

.. math::

    \newcommand{\tmp}[1]{
        J
        \vel{#1}
        \left\{
            \frac{1}{J}
            \dif{
                \left(
                    \ave{
                        \frac{J}{\sfact{#1}} \vel{1}
                    }{\gcs{#1}}
                    \ave{\vel{#1}}{\gcs{1}}
                \right)
            }{\gcs{1}}
            +
            \frac{1}{J}
            \dif{
                \left(
                    \ave{
                        \frac{J}{\sfact{#1}} \vel{2}
                    }{\gcs{#1}}
                    \ave{\vel{#1}}{\gcs{2}}
                \right)
            }{\gcs{2}}
            +
            \frac{1}{J}
            \dif{
                \left(
                    \ave{
                        \frac{J}{\sfact{#1}} \vel{3}
                    }{\gcs{#1}}
                    \ave{\vel{#1}}{\gcs{3}}
                \right)
            }{\gcs{3}}
        \right\}
    }
    \sumzc
    \sumyc
    \sumxf
    \tmp{1},

    \sumzc
    \sumyf
    \sumxc
    \tmp{2},

    \sumzf
    \sumyc
    \sumxc
    \tmp{3}.

The resulting terms are product of the velocity and the advective terms of the divergence form whose sign is flipped, and as a result the aforementioned conclusion is obtained.

*******************
Squared Temperature
*******************

We consider

.. math::

    \dtempadv{1}
    \dtempadv{2}
    \dtempadv{3}
    -
    \frac{1}{J}
    T
    \left\{
        \dif{
            \left(
                \frac{J}{\sfact{1}}
                \vel{1}
            \right)
        }{\gcs{1}}
        +
        \dif{
            \left(
                \frac{J}{\sfact{2}}
                \vel{2}
            \right)
        }{\gcs{2}}
        +
        \dif{
            \left(
                \frac{J}{\sfact{3}}
                \vel{3}
            \right)
        }{\gcs{3}}
    \right\},

giving

.. math::

    -
    \frac{1}{J}
    \dif{}{\gcs{1}}
    \left(
        \frac{J}{\sfact{1}}
        \vel{1}
        \ave{
            T
        }{\gcs{1}}
    \right)
    -
    \frac{1}{J}
    \dif{}{\gcs{2}}
    \left(
        \frac{J}{\sfact{2}}
        \vel{2}
        \ave{
            T
        }{\gcs{2}}
    \right)
    -
    \frac{1}{J}
    \dif{}{\gcs{3}}
    \left(
        \frac{J}{\sfact{3}}
        \vel{3}
        \ave{
            T
        }{\gcs{3}}
    \right).

Now we focus on the global balance of the quadratic quantity:

.. math::

    \sumzc
    \sumyc
    \sumxc
    J
    T
    \left(
        \dtempadv{1}
        \dtempadv{2}
        \dtempadv{3}
    \right),

giving

.. math::

    \newcommand{\tmp}[1]{
        \frac{1}{J}
        \dif{}{\gcs{#1}}
        \left(
            \frac{J}{\sfact{#1}}
            \vel{#1}
            \ave{T}{\gcs{#1}}
        \right)
    }
    \sumzc
    \sumyc
    \sumxc
    J
    T
    \left\{
        \tmp{1}
        +
        \tmp{2}
        +
        \tmp{3}
    \right\}.

Thus the advective contribution vanishes.

