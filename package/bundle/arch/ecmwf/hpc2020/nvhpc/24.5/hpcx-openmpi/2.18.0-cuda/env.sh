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
module_load prgenv/nvidia
module_load nvidia/24.5
module_load hpcx-openmpi/2.18.0-cuda
module_load openblas/0.3.26
# Don't load these modules if env.sh is used as part of the IFS runtime environment - only the modules above are required
if [[ ${IFS_RUNTIME_ENV:-unset} == "unset" ]]; then
  module_load python3/3.8.8-01
  module_load fftw/3.3.10
  module_load netcdf4/4.9.2
  module_load hdf5/1.14.3
  module_load cmake/3.25.2
  module_load ninja/1.10.0
  module_load fcm/2019.05.0
  module_load aec/1.0.4
fi

# MKL envs
export MKL_CBWR=AUTO,STRICT
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE # Using capital letters
export MKL_VERBOSE=${MKL_VERBOSE:-0} # if eq to 1, then each MKL func call as we go along will be output to ifs.out (stdout)
export KMP_DETERMINISTIC_REDUCTION=true

# Record the RPATH in the executable
export LD_RUN_PATH=$LD_LIBRARY_PATH

# Restore tracing to stored setting
{ if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null

