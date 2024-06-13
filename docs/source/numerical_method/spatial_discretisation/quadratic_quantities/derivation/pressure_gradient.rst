#######################
Pressure-gradient terms
#######################

The pressure-gradient terms:

.. math::

    -
    \frac{1}{\sfact{i}}
    \dif{p}{\gcs{i}}

contribute to the energy balance as follows:

.. math::

    \newcommand{\tmp}[1]{
        J
        \vel{#1}
        \frac{1}{\sfact{#1}}
        \dif{p}{\gcs{#1}}
        =
        \sumzc
        \sumyc
        \sumxc
        J
        p
        \frac{1}{J}
        \dif{
            \left(
                \frac{J}{\sfact{#1}}
                \vel{#1}
            \right)
        }{\gcs{#1}}
    }
    -
    \sumzc
    \sumyc
    \sumxf
    \tmp{1},

    -
    \sumzc
    \sumyf
    \sumxc
    \tmp{2},

    -
    \sumzf
    \sumyc
    \sumxc
    \tmp{3}.

The sum of these three relations is

.. math::

    \sumzc
    \sumyc
    \sumxc
    J
    p
    \left\{
        \ddiv{1}
        +
        \ddiv{2}
        +
        \ddiv{3}
    \right\},

which is zero because the component inside the wavy parentheses is :ref:`the incompressibility constraint <discrete_incompressibility>`.

