#
#
#   oifs-config.edit_me.sh
#
#
#   This script sets the environment for OpenIFS 48r1
#
#
#   Read this script using the command:
#
#   source ./oifs-config.edit_me.sh
#
#
#--- set machine specific settings -----------------------------

export OIFS_HOST="ecmwf"
export OIFS_PLATFORM="hpc2020"


#--- set principal OIFS variables ------------------------------
#
# CUSTOMIZED FOR USER NLD5163
#
# CODE  ----> stored in $PERM/openifs/oifs_code/branch_name/openifs-48r1
# DATA  ----> stored in $PERM/openifs/oifs_data/CYCLE/  
# EXPS  ----> stored in $PERM/openifs/oifs_exps/CYCLE/
#
# This file is therefore adapted to this particular workflow
#
#
#
#

branch_path="$(pwd)"         # We run this script in its directoryr
old=$IFS
IFS=$IFS"/"                  # We will split path in char "/". IFS variable is bash stuff.
array_path=()                # Initialize array to empty
for i in $branch_path; do    # We fill array with elements 
    array_path+=($i)
done   
IFS=$old 

branch_name=${array_path[4]} # In my current setting branch_name is in position 3 (from 0)
echo 
echo "Selecting the branch as: "$branch_name
export OIFS_CYCLE=48r1

#---Base code assumes openifs-48r1 and openifs-expt are installed
#---in $HOME. Either these can be changed by the user------------

export OIFS_HOME="${PERM}/openifs/oifs_code/"$branch_name"/openifs-48r1_fix_ifstests"

#---It is recommended that the openifs-expt and oifs_data dir
#---exist in a location designed for permanent storage-----------
export OIFS_EXPT="${PERM}/openifs/oifs_exps"/${OIFS_CYCLE} 
export OIFS_DATA_DIR="${HPCPERM}/openifs/oifs_data/"${OIFS_CYCLE}

#---Set the path for the arch directory. Depending on system,i.e.,
#---all libs are installed on the sytem, this is not required,
#---so set to an empty string OIFS_ARCH=""
#export OIFS_ARCH="./arch/ecmwf/hpc2020/gnu/"
export OIFS_ARCH="./arch/ecmwf/hpc2020/"

#---Path to the executable for 3d global model. This is the
#---default path for the exe, produced by openifs-test.sh.
#---DP means double precision. To run single precision change
#---DP to SP
export OIFS_EXEC="${OIFS_HOME}/build/bin/ifsMASTER.DP"

#---Default assumed paths, only change if you know what you are doing
export OIFS_TEST="${OIFS_HOME}/scripts/build_test"
export OIFS_RUN_SCRIPT="${OIFS_HOME}/scripts/exp_3d"
export OIFS_LOGFILE="${OIFS_HOME}/oifs_test_log.txt"

alias oenv="env -0 | sort -z | tr '\0' '\n' | grep -a OIFS_"

echo -e "\nOpenIFS environment variables are:"
echo "------------------------------------------------------"
env -0 | sort -z | tr '\0' '\n' | grep -a OIFS_
echo

#---Path to the executable for the SCM. This is the
#---default path for the exe, produced by openifs-test.sh.
#---DP means double precision. To run single precision change
#---DP to SP
export SCM_EXEC="${OIFS_HOME}/build/bin/MASTER_scm.DP"

#---Default assumed paths, only change if you know what you are doing
export SCM_TEST="${OIFS_HOME}/scripts/scm"
export SCM_RUNDIR="${OIFS_EXPT}/scm_openifs/48r1/scm-projects/ref48r1"
export SCM_PROJDIR="${OIFS_EXPT}/scm_openifs/48r1/scm-projects"
export SCM_VERSIONDIR="${OIFS_EXPT}/scm_openifs/48r1"
export SCM_LOGFILE="${SCM_RUNDIR}/scm_run_log.txt"

alias scm_env="env -0 | sort -z | tr '\0' '\n' | grep -a SCM_"

echo -e "\nSCM environment variables are:"
echo "------------------------------------------------------"
env -0 | sort -z | tr '\0' '\n' | grep -a SCM_
echo
