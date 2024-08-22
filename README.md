# Simple Navier-Stokes Solver

[![License](https://img.shields.io/github/license/NaokiHori/SimpleNSSolver)](https://opensource.org/license/MIT) 
[![LastCommit](https://img.shields.io/github/last-commit/NaokiHori/SimpleNSSolver/main)](https://github.com/NaokiHori/SimpleNSSolver/commits/main)

[![CI](https://github.com/NaokiHori/SimpleNSSolver/actions/workflows/ci.yml/badge.svg)](https://github.com/NaokiHori/SimpleNSSolver/actions/workflows/ci.yml)
[![Doc](https://github.com/NaokiHori/SimpleNSSolver/actions/workflows/documentation.yml/badge.svg)](https://naokihori.github.io/SimpleNSSolver)

[![Thumbnail](https://github.com/NaokiHori/SimpleNSSolver/blob/main/docs/source/thumbnail.gif)](https://youtu.be/WUfq8PcEhpU)

![Thumbnail](https://github.com/NaokiHori/SimpleNSSolver/blob/main/docs/source/snapshot3d.png)

## Overview

This library numerically solves the incompressible Navier-Stokes equations (coupled with a temperature field) in two- and three-dimensional Cartesian domains using the finite-difference method.

The main objective is to develop a library where the implementation and the background knowledge are closely linked via a [documentation](https://naokihori.github.io/SimpleNSSolver) and various [examples](https://naokihori.github.io/SimpleNSSolver/examples/main.html), so that users can understand *how* and *why* things are treated.

## Features

- An energy-consistent treatment of advective, pressure-gradient, and diffusive terms, correctly replicating properties of the conservation laws.
- [MPI parallelisation](https://github.com/NaokiHori/SimpleDecomp).
- Efficient FFT-based direct Poisson solver.
- Explicit / implicit treatments of diffusive terms in all spatial directions.

Please refer to the [documentation](https://naokihori.github.io/SimpleNSSolver) for details.

## Dependency

- [C compiler](https://gcc.gnu.org)
- [GNU Make](https://www.gnu.org/software/make/)
- [MPI](https://www.open-mpi.org)
- [FFTW3](https://www.fftw.org)
- [Git](https://git-scm.com)
- [Python3](https://www.python.org) with [NumPy](https://numpy.org) (for flow-field initialisation and for post-processing)

### Ubuntu

It should be convenient to use a proper package manager, e.g.:

```bash
sudo apt-get -y update
sudo apt-get -y install gcc libopenmpi-dev libfftw3-dev make
```

Also install [Python3](https://www.python.org/downloads/).

### MacOS

Installation of the `Command Line Tools for Xcode` is usually required, which is followed by

```bash
brew install gcc open-mpi fftw make
```

Also install [Python3](https://www.python.org/downloads/).

### Windows

Not supported.
Please consider to use [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/) for instance.

## Quick start

### Pre-processing

1. Prepare workplace

    ```bash
    mkdir -p /path/to/your/directory
    cd       /path/to/your/directory
    ```

1. Get source

    ```bash
    git clone --recurse-submodules https://github.com/NaokiHori/SimpleNSSolver
    cd SimpleNSSolver
    ```

1. Set initial condition

    Here Python3 is used to initialise the flow fields conveniently.
    One can give `NPY` files in different way under `initial_condition/output/`.

    ```bash
    cd initial_condition
    make output
    bash exec.sh
    cd ..
    ```

1. Build NS solver

    ```bash
    make output
    make all
    ```

### Main

```bash
bash exec.sh
```

launches the simulator and integrate the equations in time, giving e.g.

```
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
```

You see that the solver (e.g. `DOMAIN` and `FLUID`) is initialised and parameters are loaded from the `NPY` files prepared in the previous step, which is followed by the integration of the equations in time.

### Post-processing

Several log files, snapshots of the flow fields (which are used to restart the simulation and to process the flow fields later), and collected statistics are stored in `output` directory:

```
output
├── log
│  ├── xxxxx.dat
│  ├── yyyyy.dat
...
│  └── zzzzz.dat
├── save
│  ├── step00000xxxxx
│  ├── step00000yyyyy
...
│  └── step00000zzzzz
└── stat
  └── step00000zzzzz
```

Log files (files under `output/log` directory) are written in ASCII format, which are to monitor the progress.

For example, since I adopt the FFT-based Poisson solver in this project, local divergence of the flow field should be small enough, which is written in `output/log/max_divergence.dat`:

![local-divergence](https://raw.githubusercontent.com/NaokiHori/SimpleNSSolver/artifacts/artifacts/typical-2d/divergence.png)

Energy injections and dissipations are also monitored, from which the Nusselt number (computed based on several different definitions) can be extracted:

![energy-balance](https://raw.githubusercontent.com/NaokiHori/SimpleNSSolver/artifacts/artifacts/typical-2d/nusselt_time.png)

Flow fields and statistical data are stored in [NPY format](https://numpy.org/doc/stable/reference/generated/numpy.lib.format.html) using [SimpleNpyIO](https://github.com/NaokiHori/SimpleNpyIO).
When `Python3` with `NumPy` and `Matplotlib` is installed, one can easily visualise the flow fields:

![snapshot-2d](https://raw.githubusercontent.com/NaokiHori/SimpleNSSolver/artifacts/artifacts/typical-2d/snapshot.png)

or statistics (e.g., mean advective and diffusive heat transfer):

![snapshot-3d](https://raw.githubusercontent.com/NaokiHori/SimpleNSSolver/artifacts/artifacts/typical-2d/nusselt_x.png)

Note that all the results shown here are automatically updated to maintain / improve the code quality, and all scripts to produce the above figures are available in the [examples](https://naokihori.github.io/SimpleNSSolver/examples/main.html).
See the [documentation](https://naokihori.github.io/SimpleNSSolver) for more details.

## 3D simulation

By default, this project simulates two-dimensional cases because they are easy to test and thus can be a good starting point.
When a three-dimensional version is needed, checkout `3d` branch.
Note that the `main` branch contains both dimensions, which is to maintain both cases at the same time (mainly for personal use).

Please refer to the [examples](), where several small-scale 3D simulations are attempted as a part of the continuous integration.

![snapshot-3d](https://raw.githubusercontent.com/NaokiHori/SimpleNSSolver/artifacts/artifacts/typical-3d/snapshot.png)

## Contributing

Feel free to ask questions, report bugs, suggest new features, polish documentation at [issues](https://github.com/NaokiHori/SimpleNSSolver/issues).

## Further simper solver

Despite its title, this library is not as simple as it may seem, primarily due to process parallelisation and its support for both two-dimensional and three-dimensional domains.
Additionally, adding new features to this solver is challenging due to the many factors that must be taken into account, such as grid size varying from process to process, implicit time marcher, Runge-Kutta method.
To address these challenges, I have developed [another library](https://github.com/NaokiHori/NS-Sandbox) that is better suited for quickly testing new ideas and features.

## Acknowledgements

The development of this CFD solver is largely motivated by [CaNS](https://github.com/CaNS-World/CaNS) and [AFiD](https://stevensrjam.github.io/Website/afid.html).

I would like to thank [Dr. Pedro Costa](https://p-costa.github.io), [Dr. Marco Rosti](https://groups.oist.jp/cffu/marco-edoardo-rosti) and [Dr. Chris Howland](https://chowland.github.io), among others, for fruitful discussions during my time at `KTH Royal Institute of Technology in Stockholm`, `the University of Tokyo`, and `University of Twente`.

