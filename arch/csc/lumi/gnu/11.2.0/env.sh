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
module_load LUMI/22.12
module_load partition/C
module_load EasyBuild-user

module_load PrgEnv-gnu/8.3.3
module_load gcc/11.2.0
module_load cray-mpich/8.1.23
module_load cray-libsci

module_load cray-fftw/3.3.10.3
module_load cray-hdf5/1.12.2.1
module_load cray-netcdf/4.9.0.1
#module_load eigen/3.3.7
#module_load cmake/3.20.2
#module_load ninja/1.10.0
#module_load fcm/2019.05.0
module_load libaec/1.0.6-cpeGNU-22.12
  

# Setting required for bit reproducibility with Intel MKL:
export MKL_CBWR=AUTO,STRICT

# Record the RPATH in the executable
export LD_RUN_PATH=$LD_LIBRARY_PATH

# Restore tracing to stored setting
{ if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null
