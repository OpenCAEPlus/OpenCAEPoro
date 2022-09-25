# OpenCAEPoro

[![Build and publish doc](https://github.com/FaspDevTeam/OpenCAEPoro/actions/workflows/doc.yml/badge.svg)](https://github.com/FaspDevTeam/OpenCAEPoro/actions/workflows/doc.yml)
[![Github pages](https://github.com/FaspDevTeam/OpenCAEPoro/actions/workflows/gh-page.yml/badge.svg)](https://github.com/FaspDevTeam/OpenCAEPoro/actions/workflows/gh-page.yml)




| Build | Test |
|:-----:|:----:|
|[![Linux GNU Build](https://github.com/FaspDevTeam/OpenCAEPoro/actions/workflows/linux_gnu_build.yml/badge.svg)](https://github.com/FaspDevTeam/OpenCAEPoro/actions/workflows/linux_gnu_build.yml)|      |


OpenCAEPoro or OCP is part of the OpenCAEPlus project written in C++. OCP
focuses on simulating multicomponent multiphase flows in porous media. For 
more information, please see 
[OCP repository](https://faspdevteam.github.io/OpenCAEPoro/).

## Get Started

Check user manual and API's (class references from Doxygen) [website](https://porous.opencaeplus.org).

For full documentation generated by Doyxgen, see [here](https://faspdevteam.github.io/OpenCAEPoro/).

## Install
There is a top level cmake configuration file to build the OpenCAEPoro lib 
and the associate test programs suite. You can use a cmake-style approach to 
compile the package; see [the official webpage](https://cmake.org) on how to
use cmake for your own operating system. 

Before building OpenCAEPoro, you need to make sure that BLAS, LAPACK, and FASP
are available. BLAS and LAPACK are ready on most systems. The FASP package (only 
the open-source part [faspsolver](https://github.com/FaspDevTeam/faspsolver) is 
required) can be downloaded from its GitHub repository. More solver options are
included in [fasp4blkoil](https://github.com/FaspDevTeam/fasp4blkoil).

The typical command for compiling OpenCAEPoro is:

Config the environment for building with cmake:
```bash
  >>> mkdir Build; cd Build; cmake ..
```

After successfully configuring, to make the library as well as examples without
installing them, just run:
```bash
  >>> make
```

To make the library and install it, run:
```bash
  >>> make install
```

Standard **uninstall** and **clean** targets are also provided. You may safely 
remove the **Build** directory as well. 

You may also use the provided scripts to build the whole project.
```bash
  >>> source ./scripts/build/linux/gnu-build.sh Debug # change the arguement to Release for release build
  >>> source ./scripts/build/linux/clean.sh gnu Debug # For clean up the build directory
```

Or if you have pnpm installed, you can also use the package.json scripts.
```bash
  >>> pnpm build:linux:gnu:build
```

## Structure
The directory structure of OpenCAEPoro is designed as follows:
  - data/: Output files for comparison purposes
  - doc/: Documentation website
  - examples/: Input files for test examples 
  - include/: Header files
  - scripts/: Automation scripts to make life easier
  - src/: Source files
  - stash/: Files that are no longer needed, but keep for future references
  - main/: Main source code for executables 
  - CMakeLists.txt: Main cmake script
  - CONTRIBUTE.md: Guidance for open source contributors
  - LICENSE: License agreement
  - README.md: This document
  - .npmrc, pnpm-lock.yaml, pnpm-workspace.yaml, package.json: Files for the OCP website
  - .clang-format: For automatic source code formatting

## License
This software is free software distributed under the Lesser General Public
License or LGPL, version 3.0 or any later versions. This software distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with OpenCAEPoro. If not, see <http://www.gnu.org/licenses/>.
