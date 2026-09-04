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

export USER_SPACK_ROOT=/fmi/projappl/project_2003011/spack-intel
module_load spack/v0.18-user

module_load intel-oneapi-compilers-classic/2021.6.0
#module_load intel-oneapi-compilers/2022.1.0
module_load intel-oneapi-tbb/2021.6.0
module_load intel-oneapi-mpi/2021.6.0
module_load intel-oneapi-mkl/2022.1.0

module_load fftw/3.3.10-mpi
module_load netcdf-c/4.8.1
module_load netcdf-fortran/4.5.4
module_load hdf5/1.12.2-mpi
module_load cmake/3.23.1

#module_load geoconda
source /fmi/projappl/project_2003011/bergmant/openifs-48r1-pls/.oifspy/bin/activate
lspack="/appl/spack/v018/install-tree/intel-2021.6.0"
llibaec=$lspack/libaec-1.0.6-txsq2w

for lib in $llibaec; do
    export CPATH=$lib/include:$CPATH
    export LIBRARY_PATH=$lib/lib64:$LIBRARY_PATH
    export PATH=$lib/bin:$PATH
    export PKG_CONFIG_PATH=$lib/lib/pkgconfig:$PKG_CONFIG_PATH
    export CMAKE_PREFIX_PATH=$lib/.:$CMAKE_PREFIX_PATH
    export LD_LIBRARY_PATH=$lib/lib64:$LD_LIBRARY_PATH
    export LIBAEC_INSTALL_ROOT=$lib
done

# Setting required for bit reproducibility with Intel MKL:
export MKL_CBWR=AUTO,STRICT

# Record the RPATH in the executable
export LD_RUN_PATH=$LD_LIBRARY_PATH

# Restore tracing to stored setting
{ if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null

