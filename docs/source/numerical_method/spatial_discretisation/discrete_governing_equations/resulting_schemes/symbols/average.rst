#######
Average
#######

At :math:`\gcs{1}` cell faces:

.. math::

    \vat{
        \ave{q}{\gcs{1}}
    }{
        i + \frac{1}{2}
    }
    =
    \left\{
        \begin{alignedat}{2}
            & \text{Negative wall:} & \vat{q}{\frac{1}{2}}, \\
            & \text{Positive wall:} & \vat{q}{\ngp{1} + \frac{1}{2}}, \\
            & \text{Otherwise:}     & \frac{1}{2} \vat{q}{i} + \frac{1}{2} \vat{q}{i + 1}.
        \end{alignedat}
    \right.

At :math:`\gcs{1}` cell centers:

.. math::

    \vat{
        \ave{q}{\gcs{1}}
    }{
        i
    }
    =
    \frac{1}{2} \vat{q}{i - \frac{1}{2}}
    +
    \frac{1}{2} \vat{q}{i + \frac{1}{2}}.

At :math:`\gcs{2}` cell faces:

.. math::

    \vat{
        \ave{q}{\gcs{2}}
    }{
        j + \frac{1}{2}
    }
    =
    \frac{1}{2} \vat{q}{j}
    +
    \frac{1}{2} \vat{q}{j + 1}.

At :math:`\gcs{2}` cell centers:

.. math::

    \vat{
        \ave{q}{\gcs{2}}
    }{
        j
    }
    =
    \frac{1}{2} \vat{q}{j - \frac{1}{2}}
    +
    \frac{1}{2} \vat{q}{j + \frac{1}{2}}.

At :math:`\gcs{3}` cell faces:

.. math::

    \vat{
        \ave{q}{\gcs{3}}
    }{
        k + \frac{1}{2}
    }
    =
    \frac{1}{2} \vat{q}{k}
    +
    \frac{1}{2} \vat{q}{k + 1}.

At :math:`\gcs{3}` cell centers:

.. math::

    \vat{
        \ave{q}{\gcs{3}}
    }{
        k
    }
    =
    \frac{1}{2} \vat{q}{k - \frac{1}{2}}
    +
    \frac{1}{2} \vat{q}{k + \frac{1}{2}}.

