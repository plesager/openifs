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
module_purge

# Load modules
module_load gcc/11.3.0

module_load openmpi/4.1.4
module_load intel-oneapi-mkl/2022.1.0
module_load fftw/3.3.10-mpi
module_load netcdf-fortran/4.5.4
module_load netcdf-c/4.8.1
module_load hdf5/1.10.8
module_load cmake/3.23.1
module_load libaec/1.0.6

# Correct python version and libraries
module_load geoconda  

# Setting required for bit reproducibility with Intel MKL:
export MKL_CBWR=AUTO,STRICT

# Record the RPATH in the executable
export LD_RUN_PATH=$LD_LIBRARY_PATH

# Restore tracing to stored setting
{ if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null
