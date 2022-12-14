# OpenCAEPoro

[![Build and publish doc](https://github.com/OpenCAEPlus/OpenCAEPoro/actions/workflows/doc.yml/badge.svg)](https://github.com/OpenCAEPlus/OpenCAEPoro/actions/workflows/doc.yml)
[![Github pages](https://github.com/OpenCAEPlus/OpenCAEPoro/actions/workflows/gh-page.yml/badge.svg)](https://github.com/OpenCAEPlus/OpenCAEPoro/actions/workflows/gh-page.yml)

| Build | Test |
|:------|:-----|
|[![Linux GNU Build](https://github.com/OpenCAEPlus/OpenCAEPoro/actions/workflows/linux_gnu_build.yml/badge.svg)](https://github.com/OpenCAEPlus/OpenCAEPoro/actions/workflows/linux_gnu_build.yml)| [![Linux GNU Test](https://github.com/OpenCAEPlus/OpenCAEPoro/actions/workflows/linux_gnu_test.yml/badge.svg)](https://github.com/OpenCAEPlus/OpenCAEPoro/actions/workflows/linux_gnu_test.yml) |
|[![Linux Intel Build](https://github.com/OpenCAEPlus/OpenCAEPoro/actions/workflows/linux_intel_build.yml/badge.svg)](https://github.com/OpenCAEPlus/OpenCAEPoro/actions/workflows/linux_intel_build.yml)|[![Linux Intel Test](https://github.com/OpenCAEPlus/OpenCAEPoro/actions/workflows/linux_intel_test.yml/badge.svg)](https://github.com/OpenCAEPlus/OpenCAEPoro/actions/workflows/linux_intel_test.yml)|
|[![Win64 Intel Build](https://github.com/OpenCAEPlus/OpenCAEPoro/actions/workflows/windows_intel_build.yml/badge.svg)](https://github.com/OpenCAEPlus/OpenCAEPoro/actions/workflows/windows_intel_build.yml)| [![Win64 Intel Test](https://github.com/OpenCAEPlus/OpenCAEPoro/actions/workflows/windows_intel_test.yml/badge.svg)](https://github.com/OpenCAEPlus/OpenCAEPoro/actions/workflows/windows_intel_test.yml) |
|[![Mac64 Intel Build](https://github.com/OpenCAEPlus/OpenCAEPoro/actions/workflows/mac_x64_intel_build.yml/badge.svg)](https://github.com/OpenCAEPlus/OpenCAEPoro/actions/workflows/mac_x64_intel_build.yml)| [![Mac64 Intel Test](https://github.com/OpenCAEPlus/OpenCAEPoro/actions/workflows/mac_x64_intel_test.yml/badge.svg)](https://github.com/OpenCAEPlus/OpenCAEPoro/actions/workflows/mac_x64_intel_test.yml) |

## Overview

OpenCAEPoro or OCP is part of the [OpenCAEPlus](https://opencaeplus.org/) project. OCP, written in C++, focuses on simulating multicomponent multiphase flows in porous media. For more information, please refer to the [website](https://porous.opencaxplus.org).

## Get Started

For user manual and API's (class references from Doxygen), see [here](https://porous.opencaeplus.org).

For full documentation generated by Doxygen, see [here](https://opencaeplus.github.io/OpenCAEPoro/).

## Quick Installation Guide

There is a top level cmake configuration file to build the OpenCAEPoro lib and the associate test program suite. You can use a cmake-style approach to compile the package; see [the official webpage](https://cmake.org) on how to use cmake for your own operating system and tool chain.

> Note: Before building OpenCAEPoro, you need to make sure that BLAS, LAPACK, and FASP are available. BLAS and LAPACK are ready on most systems and can be found by cmake automatically. The FASP package (open-source) can be downloaded from its GitHub repository [faspsolver](https://github.com/FaspDevTeam/faspsolver).

The typical two-step `cmake` building is adopted by OpenCAEPoro:

**Step 1.** To config the environment for building with cmake:

```bash
  > mkdir Build; cd Build; cmake ..
```

**Step 2.** To build the library and test program and then install:

```bash
  > make install
```

If you do not want to `install` the lib and test examples, just run `make` without any specific target. Standard `uninstall` and `clean` targets are provided in the generated Makefile. You may safely remove the `Build` directory as well.

> Note: You can also use the provided scripts to build the whole project:

```bash
  > chmod 755 cli
  > ./cli build -c intel -t all -b Debug 
```

> Note: You may change the arguments for your own setting. For example, change `Debug` to `Release` for a release build.

## Optional Dependencies

For performance and visualization purposes, you might want to link OCP with a few optional packages, like direct solvers, multistage preconditioners, vtk, etc. One or multiple of these packages can be linked with OCP as long as they do not conflict with each other.

The required and optional external dependencies for the OpenCAEPoro library are handled by the CMakeLists file in the `external` folder. For each dependency, there is a corresponding cmake include file in the `modules` folder, dealing with things such as fetching content, creating imported targets, linking to OCP, etc.

We now give a list of optional dependencies:

### MUMPS

[MUMPS](http://mumps-solver.org/) is a massively parallel sparse direct solver mainly based on the multifrontal approach. In order to be able to call the solvers from MUMPS, you need to run:

```bash
  > mkdir Build; cd Build
  > cmake -DUSE_MUMPS=ON ..
  > make -j 8 install
```

### SuperLU

[SuperLU](https://portal.nersc.gov/project/sparse/superlu/) is a general purpose library for the LU and ILU factorization methods for linear algebraic equations. It uses MPI, OpenMP, and CUDA to support various forms of parallelism. In order to be able to call the solvers from SuperLU, you need to run:

```bash
  > mkdir Build; cd Build
  > cmake -DUSE_SUPERLU=ON ..
  > make -j 8 install
```

### Strumpack

[Strumpack](https://github.com/pghysels/STRUMPACK) is a software library providing linear algebra routines and linear system solvers for sparse and for dense rank-structured linear systems. In order to be able to call the solvers from Strumpack, you need to run:

```bash
  > mkdir Build; cd Build; 
  > cmake -DUSE_STRUMPACK=ON ..
  > make -j 8 install
```

### UMFPack in SuiteSparse

UMFPack is one of the sparse direct solvers provided in [SuiteSparse](http://suitesparse.com). In order to be able to call the solvers from UMFPack, you need to run:

```bash
  > mkdir Build; cd Build
  > cmake -DUSE_UMFPACK=ON ..
  > make -j 8 install
```

### Intel MKL Pardiso

Pardiso has two versions. We now use the version in [Intel oneMKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html#base-kit) for simplicity and ease to use with Intel compiler tool chain. In order to link with Paridso, make sure that oneAPI base toolkit has been installed and the appropriate environment variables has been set. You need to run:

```bash
    > mkdir Build; cd Build; 
    > cmake -DUSE_PARDISO=ON ..
    > make -j 8 install
```

### fasp4blkoil

The above methods are all based on Gaussian elimination. If you need efficient iterative solvers, you can link with [fasp4blkoil](https://github.com/FaspDevTeam/fasp4blkoil) or [faspcpr](https://github.com/FaspDevTeam/faspcpr). These packages contain multistage preconditioners. If you wish to use preconditioners in fasp4blkoil (with UMFPACK support), for example, you need to run:

```bash
  > mkdir Build; cd Build
  > cmake -DUSE_FASP4BLKOIL=ON -DUSE_UMFPACK=ON ..
  > make -j 8 install
```

> Note: Other direct solvers listed above can also be used with `fasp4blkoil` or `faspcpr`. Just make sure you compile them with the corresponding direct solvers.

### VTK

[VTK](https://vtk.org/) is an open-source, freely available software system for 3D computer graphics, modeling, image processing, volume rendering, scientific visualization, and 2D plotting. If you need it for visualization, you can run:

```bash
  > mkdir Build; cd Build
  > cmake -DUSE_VTK=ON ..
  > make -j 8 install
```

> Note: VTK is a very large library which may take a long time to build. Thus, though the build system will try to download and build vtk if it is not found in your system, we highly recommend you install and set up VTK on your system before building OpenCAEPoro.

### ECL

[ECL](https://github.com/equinor/ecl) a package for reading and writing the result files from the Eclipse reservoir simulator, which is often used as industry benchmark. In order to output simulation results in the Eclipse format, you need to run:

```bash
  > mkdir Build; cd Build
  > cmake -DUSE_ECL=ON ..
  > make -j 8 install
```

## Structure

The directory structure of OpenCAEPoro is designed as follows:

- `data/` : Output files for comparison purposes
- `doc/` : Documentation website
- `examples/` : Input files for test examples
- `external/` : External dependencies
- `include/` : Header files
- `regression/` : Regression tests
- `src/` : Source files
- `stash/` : Files that are no longer needed, keep for future references
- `main/` : Main source code for executables
- `modules/` : Cmake files for finding and setting dependencies
- CMakeLists.txt: Main cmake script
- CMakePresets.json: Preset settings for cmake
- CONTRIBUTE.md: Guidance for open source contributors
- LICENSE: License agreement
- README.md: This document
- cli, cli.bat: Command line interface for build, test, and doc
- .clang-format: For automatic source code formatting

## Test

There are three levels of tests:

- Unit test: for testing individual code components
- Integrate test: for testing interfaces between components
- System test: for testing and benchmark of the whole program  

Currently, only system tests are provided in the `examples/` folder.

## License

This software is free software distributed under the Lesser General Public
License or LGPL, version 3.0 or any later versions. This software distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with OpenCAEPoro. If not, see <http://www.gnu.org/licenses/>.
