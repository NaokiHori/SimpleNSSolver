#################
initial_condition
#################

********
Overview
********

This directory contains a Python script to initialise the domain and the flow field, which will be simulated by the main solver.
Note that this initialiser is not parallelised for simplicity.

*************
Configuration
*************

See ``main.py``.

*****
Usage
*****

.. code-block:: console

   make output
   bash exec.sh

giving several ``NPY`` files under the specified directory (``output`` by default).

These files will be loaded by the main simulator.

