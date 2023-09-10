
.. include:: /references.txt

######################
Grossmann-Lohse theory
######################

In this section, I investigate the relation between :math:`Ra` and :math:`Nu`, which is known as the Grossmann-Lohse theory (|GROSSMANN2000|).

*************
Configuration
*************

The Prandtl number is fixed to :math:`0.7`, and the Rayleigh number is varied from :math:`10^4` to :math:`10^8` to see the relation between :math:`Ra` and :math:`Nu`.

.. mydetails:: Details

   .. literalinclude:: data/exec.sh
      :language: sh

******
Result
******

.. image:: data/nu_ra.png
   :width: 600

The black-dashed line shows

.. math::

   Nu \propto Ra^{0.28},

which is reported by e.g. |KOOLOTH2021|.

.. mydetails:: Script

   .. literalinclude:: data/process.py
      :language: python
      :linenos:

