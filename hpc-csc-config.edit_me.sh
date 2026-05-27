#
#
#   hpc-config.edit_me.sh
#
#
#   This script sets the HPC environment for CSC platforms (puhti,
#   lumi, and mahti) for OpenIFS 48r1
#
#   The script is ran in oifs-config.edit_me.sh, which requires
#   editing for that purpose.
#
#

#--- set default values to unset for some variables ----------------------
unset HPC_PARTITION
unset HPC_ACCOUNT
unset DEFAULT_NUM_THREADS
unset IGT_BUILD_LAUNCHER
unset IGT_TEST_LAUNCHER
unset HPC_FLAGS

#--- set HPC architecture specific settings -----------------------------

if [[ $HPC_PLATFORM == "hpc2020" ]]; then
    #--- HPC2020 SPECIFIC SETTINGS

    #--- global variables
    export DEFAULT_NUM_THREADS=64

    #---
    
elif [[ $HPC_PLATFORM == "puhti" ]]; then
    #--- PUHTI SPECIFIC SETTINGS
    #--- local variables
    PARTITION="--partition=fmi"
    ACCOUNT="--account=project_2003011"
    MEM="--mem=0"

    # Add exclusive flag for Puhti:
    # start a single tasks reserving a full node (shared otherwise),
    # also change the default TMPDIR (changed in openifs-bundle)
    HPC_FLAGS="-n 1 --exclusive --export=ALL,MY_TMP=/dev/shm"

    #--- global variables
    # Since we need to reserve the full node, make use of it
    export DEFAULT_NUM_THREADS=40

    # Make sure GIT is loaded
    module load git

    # Overwrite compile flags. There's a bug in Puhti installed libxml
    # resulting in floating point exception due to the use of -ffpe-trap=
    cp -f /fmi/projappl/project_2003011/bergmant/openifs-48r1-pls/arch/csc/puhti/compile_flags_puhti.cmake /fmi/projappl/project_2003011/bergmant/openifs-48r1-pls/ifs-source/cmake/compile_flags.cmake

    #---
    
elif [[ $HPC_PLATFORM == "lumi" ]]; then
    #--- LUMI SPECIFIC SETTINGS
    #--- local variables
    PARTITION="--partition=small"
    ACCOUNT="--account=project_462000178"
    MEM="--mem=60GB"

    #--- global variables
    export DEFAULT_NUM_THREADS=16
    export EBU_USER_PREFIX=/project/project_462000178/EasyBuild

    #---

elif [[ $HPC_PLATFORM == "mahti" ]]; then
    #--- MAHTI SPECIFIC SETTINGS
    #--- local variables
    PARTITION="--partition=medium"
    ACCOUNT="--account=project_2001029"
    MEM="--mem=60GB"

    # Need specify the number of tasks for slurm,
    # otherwise thread number used as tasks
    HPC_FLAGS="-n 1"

    #--- global variables
    export DEFAULT_NUM_THREADS=10

    #---
    
else
    #--- UNSPECIFIED PLATFORM DEFAULTS
    # Note: Default threads = 8 works on Mac M1 but is
    # potentially too high for older systems
    export DEFAULT_NUM_THREADS=8
fi

#--- set non-default launcher options
if [[ $HPC_HOST == "csc" ]]; then
    export IGT_BUILD_LAUNCHER="srun $HPC_FLAGS -c ${DEFAULT_NUM_THREADS} ${MEM} --time=60 $ACCOUNT $PARTITION"
    export IGT_TEST_LAUNCHER="salloc -n 8 --mem=20GB --time=60 $ACCOUNT $PARTITION"
fi
