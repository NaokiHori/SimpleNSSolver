############
Introduction
############

********
Overview
********

This library is developed with the following three main reasons in mind.

#. Bridging the gap between methods and implementations

   Despite the existence of numerous elegant libraries, not all of them have comprehensive documentation.
   In some cases, there can be a significant discrepancy between the methods (equations) and the actual implementations (code).
   The goal of this project is to reduce this gap and provide clear links.

#. Ensuring validation and verification

   While code validations are typically shown in publications, it can be challenging to determine if the code really works, and sometimes compiling the code itself can be tedious.
   This project aims to validate the code and automatically publish the results when changes are made, ensuring that the library is always well-validated and verified.

#. Shedding light on small but non-trivial details

   Although the background algorithm of this solver (integration of mass, momentum, and scalar fields) is well-known and straightforward, minor details are often overlooked, and the correct approaches can be counter-intuitive.
   Adopting Rayleigh-Bénard convection as a model problem, this project aims to shed light on these small but non-trivial details, such as discrete energy conservation, calculation of dissipation rates, and Nusselt number agreements, providing insights into these important aspects.

************
Dependencies
************

======
Solver
======

* `C compiler <https://gcc.gnu.org>`_
* `GNU Make <https://www.gnu.org/software/make/>`_
* `MPI <https://www.open-mpi.org>`_
* `FFTW3 <https://www.fftw.org>`_

=======
Utility
=======

* `Git <https://git-scm.com>`_
* `Python <https://www.python.org>`_ with `NumPy <https://numpy.org>`_ (for easy flow initialisation)

*****
Usage
*****

#. Prepare workplace

   .. code-block:: console

      $ mkdir -p /path/to/your/directory
      $ cd       /path/to/your/directory

#. Get source

   .. code-block:: console

      $ git clone --recurse-submodules https://github.com/NaokiHori/SimpleNSSolver
      $ cd SimpleNSSolver

#. Set initial condition

   A set of initial conditions must be provided, as this is a PDE solver:

   .. code-block:: console

      $ cd initial_condition
      $ make output
      $ bash exec.sh
      $ cd ..

#. Build solver

   .. code-block:: console

      $ make output
      $ make all

#. Execute

   .. code-block:: console

      $ bash exec.sh

Run-time parameters are defined in `exec.sh <https://github.com/NaokiHori/SimpleNSSolver/blob/main/exec.sh>`_ as environment variables.

