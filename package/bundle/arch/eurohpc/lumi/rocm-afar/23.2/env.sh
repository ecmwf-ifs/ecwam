# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Source me to get the correct configure/build/run environment

# Store tracing and disable (module is *way* too verbose)
{ tracing_=${-//[^x]/}; set +x; } 2>/dev/null

module_load() {
  echo "+ module load $1"
  module load $1
}
module_unload() {
  echo "+ module unload $1"
  module unload $1
}

# Unload to be certain
module reset

# Load modules
module_load LUMI/25.09
module_load partition/G
module_load cray-mpich/9.0.1
module_load craype-network-ofi
module_load buildtools/25.09
module_load cray-python/3.11.7

LD_LIBRARY_PATH=/pfs/lustrep4/scratch/project_465000527/nawabahm/therock-afar-23.2.1-gfx90a-7.13.0-7357b5084b/lib:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=/pfs/lustrep4/scratch/project_465000527/nawabahm/therock-afar-23.2.1-gfx90a-7.13.0-7357b5084b/lib/llvm/lib:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=/pfs/lustrep4/scratch/project_465000527/nawabahm/local/libffi-3.2.1/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH
PATH=/pfs/lustrep4/scratch/project_465000527/nawabahm/therock-afar-23.2.1-gfx90a-7.13.0-7357b5084b/lib/llvm/bin:$PATH
PATH=/pfs/lustrep4/scratch/project_465000527/nawabahm/therock-afar-23.2.1-gfx90a-7.13.0-7357b5084b/bin:$PATH
export PATH=$PATH

# Export environment variable3s
export MPI_HOME=${MPICH_DIR}

export CC=amdclang CXX=amdclang++ FC=amdflang

module list

set -x

ulimit -s unlimited
ulimit -l unlimited

# Restore tracing to stored setting
{ if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null

# export ECBUILD_TOOLCHAIN="./toolchain.cmake"
path=$BASH_SOURCE
DIR_PATH=$(dirname $path)
export ECBUILD_TOOLCHAIN=$DIR_PATH/toolchain.cmake
