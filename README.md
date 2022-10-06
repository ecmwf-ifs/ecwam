ecWAM
*****

The ECMWF wave model WAM

Introduction
============


License
=======

ecWAM is distributed under the Apache License Version 2.0.
See `LICENSE` file for details.

Installing ecWAM 
================

Supported Platforms
-------------------

- Linux
- Apple MacOS

Other UNIX-like operating systems may work too out of the box.

Requirements
------------
- Fortran and C compiler, and optionally C++ compiler
- CMake (see https://cmake.org)
- ecbuild (see https://github.com/ecmwf/ecbuild)
- fiat
- eccodes

Further optional dependencies:
- MPI Fortran libraries
- multio
- NEMO

Building ecWAM
-------------

Environment variables 

    $ export ecbuild_ROOT=<path-to-ecbuild>
    $ export MPI_HOME=<path-to-MPI>
    $ export fiat_ROOT=<path-to-fiat>
    $ export CC=<path-to-C-compiler>
    $ export FC=<path-to-Fortran-compiler>
    $ export CXX=<path-to-C++-compiler> 

You must compile ecWAM out-of-source, so create a build-directory

    $ mkdir build && cd build
 
Configuration of the build happens through standard CMake

    $ cmake ..

Extra options can be added to the `cmake` command to control the build:

 - `-DCMAKE_BUILD_TYPE=<Debug|RelWithDebInfo|Release|Bit>` default=RelWithDebInfo (typically `-O2 -g`)
 - `-DENABLE_TESTS=<ON|OFF>` 
 - `-DENABLE_MPI=<ON|OFF>` 
 - `-DENABLE_OMP=<ON|OFF>`
 - `-DCMAKE_INSTALL_PREFIX=<install-prefix>`

More options to control compilation flags, only when defaults are not sufficient

 - `-DOpenMP_Fortran_FLAGS=<flags>`
 - `-DCMAKE_Fortran_FLAGS=<fortran-flags>`
 - `-DCMAKE_C_FLAGS=<c-flags>`

Once this has finished successfully, run ``make`` and ``make install``.

Optionally, tests can be run to check succesful compilation, when the feature TESTS is enabled (`-DENABLE_TESTS=ON`, default ON)

    $ ctest

Reporting Bugs
==============

TODO

Contributing
============

TODO

