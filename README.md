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
--------------

Environment variables 

    $ export ecbuild_ROOT=<path-to-ecbuild>
    $ export MPI_HOME=<path-to-MPI>
    $ export fiat_ROOT=<path-to-fiat>
    $ export eccodes_ROOT=<path-to-eccodes>
    $ export CC=<path-to-C-compiler>
    $ export FC=<path-to-Fortran-compiler>
    $ export CXX=<path-to-C++-compiler> 

If you want to pre-download or install extra data files, run in the source-directory before CMake (re)configuration:

    $ share/ecwam/data/populate.sh

You must compile ecWAM out-of-source, so create a build-directory

    $ mkdir build && cd build
 
Configuration of the build happens through standard CMake

    $ cmake <path-to-source> 

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


Running ecWAM
=============

To run, use the scripts in the build or install directory.
Outputs will be stored in "wamrun" directory

1) Configuration options

A YAML configuration file is required. See [tests](tests) directory for examples.

2) Create the bathymetry from ETOPO1 dataset

```shell
share/ecwam/scripts/ecwam_run_create_bathymetry.sh --run-dir=$(pwd)/wamrun --config=<path-to-config.yml>
```

3) Create grid tables

```shell
share/ecwam/scripts/ecwam_run_preproc.sh --run-dir=$(pwd)/wamrun --config=<path-to-config.yml>
```

4) Create initial conditions

```shell
share/ecwam/scripts/ecwam_run_preset.sh --run-dir=$(pwd)/wamrun --config=<path-to-config.yml>
```

5) Run wave model

```shell
share/ecwam/scripts/ecwam_run_model.sh --run-dir=$(pwd)/wamrun --config=<path-to-config.yml>
```

These shell scripts set up running of wave model executables.
Anything present in a environment variable `LAUNCH`
will be prefixed to the launching of the executable.
For instance to run with 4 MPI tasks you should, depending on the MPIEXEC:
```
export LAUNCH="mpirun -np 4"
```

If `--run-dir` argument is not specified, or `RUN_DIR` environment variable is not set, then the current directory is used.
If `--config` argument is not specified, it is assumed that a file called `config.yml` is pressent in the run-dir.


Reporting Bugs
==============

TODO

Contributing
============

TODO

