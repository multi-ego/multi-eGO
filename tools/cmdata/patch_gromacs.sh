#!/bin/bash

CURRENT=$(pwd)
GROMACS_ROOT_DIR=$1

[[ $# -eq 0 ]] && [[ ! -d $GROMACS_ROOT_DIR ]] &&  echo "Please supply a valid gromacs root directory as a command-line argument" && exit 0

diff -aur ${GROMACS_ROOT_DIR}/src/gromacs/trajectoryanalysis/modules.cpp modules.cpp > patch.diff
cp cmdata.cpp cmdata.h ${GROMACS_ROOT_DIR}/src/gromacs/trajectoryanalysis/modules/
touch ${GROMACS_ROOT_DIR}/src/gromacs/trajectoryanalysis/* 

patch -d/ -p0 < patch.diff
