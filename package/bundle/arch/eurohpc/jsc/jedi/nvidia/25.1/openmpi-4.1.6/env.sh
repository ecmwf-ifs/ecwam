# Source me to get the correct configure/build/run environment

# Store tracing and disable (module is *way* too verbose)
{ tracing_=${-//[^x]/}; set +x; } 2>/dev/null

module_load() {
  echo "+ module load $*"
  module load $*
}
module_unload() {
  echo "+ module unload $*"
  module unload $*
}
module_purge() {
  echo "+ module purge"
  module purge
}

# Unload all modules to be certain
[[ ${IFS_RUNTIME_ENV:-unset} == "unset" ]] && module_purge

# Load modules
module_load Stages/2024
module_load OpenSSL/1.1
module_load CUDA/12
module_load StdEnv/2024
module_load NVHPC/25.1-CUDA-12
module_load git/2.41.0-nodocs
module_load CMake/3.26.3
module_load OpenMPI/4.1.6
module_load Python/3.11.3
[[ ${ECTRANS_GPU_AWARE_MPI:-unset} == "1" ]] && module_load MPI-settings/CUDA

# Record the RPATH in the executable
export LD_RUN_PATH=$LD_LIBRARY_PATH

export FC=nvfortran
export CC=nvc
export CXX=nvc++
# Restore tracing to stored setting
{ if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null

