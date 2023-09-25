###########################
Simple Navier-Stokes Solver
###########################

|License|_ |LastCommit|_

.. |License| image:: https://img.shields.io/github/license/NaokiHori/SimpleNSSolver
.. _License: https://opensource.org/license/MIT

.. |LastCommit| image:: https://img.shields.io/github/last-commit/NaokiHori/SimpleNSSolver/main
.. _LastCommit: https://github.com/NaokiHori/SimpleNSSolver/commits/main

|CI|_ |DOC|_

.. |CI| image:: https://github.com/NaokiHori/SimpleNSSolver/actions/workflows/ci.yml/badge.svg
.. _CI: https://github.com/NaokiHori/SimpleNSSolver/actions/workflows/ci.yml

.. |DOC| image:: https://github.com/NaokiHori/SimpleNSSolver/actions/workflows/documentation.yml/badge.svg
.. _DOC: https://naokihori.github.io/SimpleNSSolver

.. thumbnails

.. image:: https://github.com/NaokiHori/SimpleNSSolver/blob/main/docs/source/thumbnail.gif
   :width: 800

.. image:: https://github.com/NaokiHori/SimpleNSSolver/blob/main/docs/source/snapshot3d.png
   :width: 800

.. shortcuts

.. _implementations: https://naokihori.github.io/SimpleNSSolver/implementation/main.html
.. _theories: https://naokihori.github.io/SimpleNSSolver/equations/main.html
.. _numerics: https://naokihori.github.io/SimpleNSSolver/numerical_method/main.html
.. _examples: https://naokihori.github.io/SimpleNSSolver/examples/main.html
.. _documentation: https://naokihori.github.io/SimpleNSSolver

.. contents::
   :depth: 1

********
Overview
********

This library numerically solves the incompressible Navier-Stokes equations (coupled with a temperature field) in two- and three-dimensional Cartesian domains using the finite-difference method.

The main objective is to develop a library where the `implementations`_ and their background knowledge (`theories`_ and `numerics`_) are closely linked via a `documentation`_ and various `examples`_, so that users can understand *how* and *why* things are treated.

********
Features
********

* Strong link between `theories`_, `numerics`_, `implementations`_ and `examples`_ through a `documentation`_.
* Automatic code validation (`examples`_) using `GitHub Actions <https://github.com/features/actions>`_.
* `Numerically perfect mass and momentum balances <https://naokihori.github.io/SimpleNSSolver/examples/typical/main.html>`_.
* `Numerically perfect Nusselt number agreements <https://naokihori.github.io/SimpleNSSolver/examples/nu/main.html>`_.
* `Conservation of quadratic quantities (squared velocity and temperature) for inviscid fluids <https://naokihori.github.io/SimpleNSSolver/examples/energy/main.html>`_ and resulting extreme stability.
* `Pencil-based MPI parallelisation <https://github.com/NaokiHori/SimpleDecomp>`_, from single process to O(10^4) processes.
* `Efficient FFT-based direct Poisson solver <https://naokihori.github.io/SimpleNSSolver/implementation/fluid/compute_potential/main.html>`_.
* `Explicit/implicit treatment of diffusive terms in all spatial directions <https://naokihori.github.io/SimpleNSSolver/implementation/linear_system.html>`_ are easily switchable.

Please refer to the `documentation`_ for details.

**********
Dependency
**********

* `C compiler <https://gcc.gnu.org>`_
* `GNU Make <https://www.gnu.org/software/make/>`_
* `MPI <https://www.open-mpi.org>`_
* `FFTW3 <https://www.fftw.org>`_
* `Git <https://git-scm.com>`_
* `Python3 <https://www.python.org>`_ with `NumPy <https://numpy.org>`_ (for flow-field initialisation and for post-processing)

======
Ubuntu
======

It should be convenient to use a proper package manager, e.g.:

.. code-block:: console

   sudo apt-get -y update
   sudo apt-get -y install gcc libopenmpi-dev libfftw3-dev make

Also install `Python3 <https://www.python.org/downloads/>`_.

=====
MacOS
=====

Installation of the ``Command Line Tools for Xcode`` is usually required, which is followed by

.. code-block:: console

   brew install gcc open-mpi fftw make

Also install `Python3 <https://www.python.org/downloads/>`_.

=======
Windows
=======

Not supported.
Please consider to use `Windows Subsystem for Linux <https://learn.microsoft.com/en-us/windows/wsl/>`_ for instance.

***********
Quick start
***********

===========
Preparation
===========

#. Prepare workplace

   .. code-block:: console

      mkdir -p /path/to/your/directory
      cd       /path/to/your/directory

#. Get source

   .. code-block:: console

      git clone --recurse-submodules https://github.com/NaokiHori/SimpleNSSolver
      cd SimpleNSSolver

#. Set initial condition

   Here ``Python3`` is used to initialise the flow fields conveniently.
   One can give ``NPY`` files in different way under ``initial_condition/output/``.

   .. code-block:: console

      cd initial_condition
      make output
      bash exec.sh
      cd ..

#. Build NS solver

   .. code-block:: console

      make output
      make all

==========
Simulation
==========

.. code-block:: console

   bash exec.sh

launches the simulator and integrate the equations in time, giving e.g.

.. code-block:: text

   DOMAIN
      glsizes[0]: 128
      glsizes[1]: 256
      lengths[0]:  1.0000000e+00
      lengths[1]:  2.0000000e+00
   FLUID
      Ra:  1.0000000e+08
      Pr:  1.0000000e+01
      Momentum    diffusivity:  3.1622777e-04
      Temperature diffusivity:  3.1622777e-05
      diffusive treatment in x: implicit
      diffusive treatment in y: explicit
   LOGGING
      next:  5.000e-01
      rate:  5.000e-01
   SAVE
      dest: output/save/step
      next:  2.000e+01
      rate:  2.000e+01
   STATISTICS
      dest: output/stat/step
      next:  1.000e+02
      rate:  1.000e-01
   step: 0, time:  0.0000000e+00
   timemax:  2.0000000e+02, wtimemax:  6.0000000e+02
   coefs: (adv)  9.500e-01, (dif)  9.500e-01
   DFT-based solver is used
   step   11, time   0.5, dt 4.58e-02, elapsed  2.1 [sec]
   step   22, time   1.0, dt 4.58e-02, elapsed  2.2 [sec]
   step   33, time   1.5, dt 4.58e-02, elapsed  2.3 [sec]
   step   44, time   2.0, dt 4.58e-02, elapsed  2.4 [sec]
   step   55, time   2.5, dt 4.58e-02, elapsed  2.4 [sec]
   ...
   step 8193, time 197.5, dt 3.06e-02, elapsed 91.9 [sec]
   step 8210, time 198.0, dt 2.79e-02, elapsed 92.2 [sec]
   step 8228, time 198.5, dt 2.79e-02, elapsed 92.5 [sec]
   step 8246, time 199.0, dt 2.90e-02, elapsed 93.0 [sec]
   step 8263, time 199.5, dt 3.07e-02, elapsed 93.2 [sec]

You see that the solver (e.g. ``DOMAIN`` and ``FLUID``) is initialised and parameters are loaded from the ``NPY`` files prepared in the previous step, which is followed by the integration of the equations in time.

===============
Post-processing
===============

Several log files, snapshots of the flow fields (which are used to restart the simulation and to process the flow fields later), and collected statistics are stored in ``output`` directory:

.. code-block:: text

   output
   ├── log
   │  ├── divergence.dat
   │  ├── energy.dat
   │  ├── momentum.dat
   │  ├── nusselt.dat
   │  └── progress.dat
   ├── save
   │  ├── step00000xxxxx
   │  ├── step00000yyyyy
   ...
   │  └── step00000zzzzz
   └── stat
      └── step00000zzzzz

Log files (files under ``output/log`` directory) are written in ASCII format, which are to monitor the progress.

For example, since I adopt the FFT-based Poisson solver in this project, local divergence of the flow field should be small enough, which is written in ``output/log/divergence.dat``:

.. image:: https://naokihori.github.io/SimpleNSSolver/_images/divergence_2d.png
   :width: 50%

Also the Nusselt numbers (computed based on several different definitions, see the `documentation`_) are monitored and written in ``output/log/nusselt.dat``:

.. image:: https://naokihori.github.io/SimpleNSSolver/_images/nusselt_time_2d.png
   :width: 50%

Flow fields and statistical data are stored in `NPY format <https://numpy.org/doc/stable/reference/generated/numpy.lib.format.html>`_ using `SimpleNpyIO <https://github.com/NaokiHori/SimpleNpyIO>`_.
When ``Python3`` with ``NumPy`` and ``Matplotlib`` is installed, one can easily visualise the flow fields:

.. image:: https://naokihori.github.io/SimpleNSSolver/_images/snapshot_2d.png
   :width: 50%

Also it is trivial to extract statistics.
For example, this plot shows the fluctuations of the flow quantities:

.. image:: https://naokihori.github.io/SimpleNSSolver/_images/std_2d.png
   :width: 50%

or the mean heat flux:

.. image:: https://naokihori.github.io/SimpleNSSolver/_images/nusselt_x_2d.png
   :width: 50%

By varying the parameter (in particular the Rayleigh number ``Ra``), one can observe a famous scaling law:

.. image:: https://naokihori.github.io/SimpleNSSolver/_images/nu_ra.png
   :width: 50%

Note that all the results shown here are automatically updated to maintain / improve the code quality, and all scripts to produce the above figures are available in the `examples`_.
See the `documentation`_ for more details.

*************
3D simulation
*************

By default, this project simulates two-dimensional cases because they are easy to test and thus can be a good starting point.
When the three-dimensional counterpart is needed, checkout the ``3d`` branch.
Note that the ``main`` branch contains both dimensions, which is for the developers to maintain both cases at the same time.

Please refer to the `examples`_, where several small-scale 3D simulations are attempted as a part of the continuous integration.

.. image:: https://naokihori.github.io/SimpleNSSolver/_images/snapshot_3d.png
   :width: 50%

************
Contributing
************

Feel free to ask questions, to report bugs, or to suggest new features at `issues <https://github.com/NaokiHori/SimpleNSSolver/issues>`_.

****************
Acknowledgements
****************

The development of this CFD solver is largely motivated by `CaNS <https://github.com/CaNS-World/CaNS>`_ and `AFiD <https://stevensrjam.github.io/Website/afid.html>`_.

I would like to thank `Dr. Pedro Costa <https://p-costa.github.io>`_, `Dr. Marco Rosti <https://groups.oist.jp/cffu/marco-edoardo-rosti>`_ and `Dr. Chris Howland <https://chowland.github.io>`_, among others, for fruitful discussions during my time at KTH Royal Institute of Technology in Stockholm, the University of Tokyo and University of Twente.

