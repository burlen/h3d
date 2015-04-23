#!/bin/bash
WORKDIR=/home/danlipsa/src/h3d-newcode/src

export DATA_DIRECTORY=$WORKDIR/data
export RESTART_DIRECTORY=$WORKDIR/restart_files
export INPUT_DIRECTORY=$WORKDIR/../inputs
export SOURCE_DIRECTORY=$WORKDIR

mkdir -p $DATA_DIRECTORY
mkdir -p $RESTART_DIRECTORY

cd $WORKDIR
cp $INPUT_DIRECTORY/finput.dat $DATA_DIRECTORY/finput.dat
#cp ${HOME}/.cleanup_status_init $DATA_DIRECTORY/.cleanup_status_init
touch $DATA_DIRECTORY/.cleanup_status_init
touch $DATA_DIRECTORY/.cleanup_status

export MPI_TYPE_MAX=65536
export MPI_REQUEST_MAX=65536

mpiexec -np 1 ./3dh > 3dhout.txt
