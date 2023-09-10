
.. _implementation:

##############
Implementation
##############

The purpose of this section is to explain the details of the implementation.

This project adopts a directory structure which is widely used:

* Under `include <https://github.com/NaokiHori/SimpleNSSolver/tree/main/include>`_, structures and non-static functions are declared.

* Under `src <https://github.com/NaokiHori/SimpleNSSolver/tree/main/src>`_, all functions are implemented.

Although I try to document in the code, it cat be too verbose if I explain all things there.
The following pages are mainly used to give additional information.
In particular I aim to describe :ref:`the discretisations of the equations <numerics>` and their implementations side-by-side.

.. toctree::
   :maxdepth: 1

   domain
   fileio
   fluid/main
   halo/main
   integrate
   linear_system
   logging/main
   param
   statistics
   tdm/main

.. seealso::

   Please check the README files placed under ``src`` to explain the role of each directory and file.

.. note::

   The main solver always try to load a flow field from a specified directory, i.e. a set of initial conditions is to be given externally.
   See `initial_condition/main.py <https://github.com/NaokiHori/SimpleNSSolver/blob/2d/initial_condition/main.py>`_.

*************
`src/main.c`_
*************

.. _src/main.c: https://github.com/NaokiHori/SimpleNSSolver/blob/main/src/main.c

This file contains a main function ``main`` which is an entry point of this solver and calls all other functions.

.. mydeclare:: /../../src/main.c
   :language: c
   :tag: main

In the following, the overview of the solver is briefly described.

=====
Setup
=====

Since this is a PDE solver, initial condition for each scalar field and the coordinates (i.e. grid points) are required, which should be stored as ``NPY`` files under a specific directory.
Normally you can create the initial condition using the separated program under directory ``initial_condition/``.
I expect the name of the directory where all the initial flow fields and the coordinate systems are stored is passed as the argument: for example,

.. code-block:: console

   mpirun -n 16 ./a.out initial_condition/output

when you start from the initial condition, while

.. code-block:: console

   mpirun -n 16 ./a.out output/save/step0000012345

when you restart from the intermediate flow field at ``12345`` step.

Run-time parameters (such as the maximum simulation time) are given as environment variables (i.e. ``export KEY=VALUE``) when the program is launched.
Although a straightforward way is to give them in front of ``mpirun`` (e.g. ``KEY1=VALUE1 KEY2=VALUE2 ... mpirun ... ./a.out``), it is cumbersome.
A more convenient way is to define them in a file: see `exec.sh <https://github.com/NaokiHori/SimpleNSSolver/blob/main/exec.sh>`_ for example.
Other parameters, which are less frequently changed but still possibly changed, are configured under :ref:`param <param>`.

==============
Initialisation
==============

First of all, ``MPI`` is launched by invoking ``MPI_Init``, and time when the simulator is started is recorded:

.. myliteralinclude:: /../../src/main.c
   :language: c
   :tag: launch MPI, start timer

``wtimes`` are used to check the elapsed time to terminate the simulator (see below).

I expect the name of the directory where the initial conditions are stored is passed as the argument:

.. myliteralinclude:: /../../src/main.c
   :language: c
   :tag: find name of directory where IC is stored

From this directory, the time step ``step`` and the time units ``time`` are loaded:

.. myliteralinclude:: /../../src/main.c
   :language: c
   :tag: initialise time step and time units

Also from this directory, the whole initial flow information is loaded and assigned to the corresponding structures:

* Domain decomposition and coordinate systems: ``domain_t``

* Flow field, velocity, pressure, temperature: ``fluid_t``

.. myliteralinclude:: /../../src/main.c
   :language: c
   :tag: initialise structures

Finally other auxiliary objects (logger, flow fields saver, statistics collector) are initialised and scheduled:

.. myliteralinclude:: /../../src/main.c
   :language: c
   :tag: initialise auxiliary objects

=========
Main loop
=========

In this part, :ref:`the main solver <fluid>` is iteratively invoked to integrate the equations in time.
Please refer to :ref:`the temporal discretisation <temporal_discretisation>` for more details about the method.

After the scalar fields have been updated, the time step ``step`` and the time units ``time`` should also be synchronised:

.. myliteralinclude:: /../../src/main.c
   :language: c
   :tag: update step and simulation / wall time

Also I check the elapsed time by recording the current wall time.

The simulator will abort if one of the abort conditions is met:

.. myliteralinclude:: /../../src/main.c
   :language: c
   :tag: terminate if one of the following conditions is met

Values to be monitored during the simulation (e.g. divergence) are computed and written to the screen or to files regularly:

.. myliteralinclude:: /../../src/main.c
   :language: c
   :tag: compute and output log regularly

Also the flow fields are regularly dumped:

.. myliteralinclude:: /../../src/main.c
   :language: c
   :tag: save flow fields regularly

In addition, temporally-averaged statistics are regularly computed and collected:

.. myliteralinclude:: /../../src/main.c
   :language: c
   :tag: collect statistics regularly

============
Finalisation
============

The last flow fields (velocity, pressure, and temperature), which are necessary to restart the simulation later, are saved:

.. myliteralinclude:: /../../src/main.c
   :language: c
   :tag: save final flow fields

Also the statistics are saved to files if collected:

.. myliteralinclude:: /../../src/main.c
   :language: c
   :tag: save collected statistics

Finally ``MPI`` is finalised:

.. myliteralinclude:: /../../src/main.c
   :language: c
   :tag: finalise MPI

