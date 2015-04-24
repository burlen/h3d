#!/bin/bash
WORKDIR=/home/utkarsh/Kitware/ParaView3/Sensei/h3d/src

export DATA_DIRECTORY=$WORKDIR/data
export RESTART_DIRECTORY=$WORKDIR/restart_files
export INPUT_DIRECTORY=$WORKDIR/../inputs
export SOURCE_DIRECTORY=$WORKDIR

mkdir -p $DATA_DIRECTORY
mkdir -p $RESTART_DIRECTORY

cd $WORKDIR
cp $INPUT_DIRECTORY/finput.small.dat $DATA_DIRECTORY/finput.dat
#cp ${HOME}/.cleanup_status_init $DATA_DIRECTORY/.cleanup_status_init
touch $DATA_DIRECTORY/.cleanup_status_init

export MPI_TYPE_MAX=65536
export MPI_REQUEST_MAX=65536

mpiexec -np 16 ./3dh > 3dhout
