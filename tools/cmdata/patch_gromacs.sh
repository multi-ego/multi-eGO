#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
GROMACS_ROOT_DIR=$1

[[ $# -eq 0 ]] && [[ ! -d $GROMACS_ROOT_DIR ]] &&  echo "Please supply a valid gromacs root directory as a command-line argument" && exit 0

diff -aur ${GROMACS_ROOT_DIR}/src/gromacs/trajectoryanalysis/modules.cpp ${SCRIPT_DIR}/modules.cpp > patch.diff
cp ${SCRIPT_DIR}/cmdata.cpp ${SCRIPT_DIR}/cmdata.h ${GROMACS_ROOT_DIR}/src/gromacs/trajectoryanalysis/modules/
touch ${GROMACS_ROOT_DIR}/src/gromacs/trajectoryanalysis/* 

patch -d/ -p0 < patch.diff
