#PBS -S /bin/csh
#PBS -N testlnewc
#PBS -o outfile
#PBS -l select=256:ncpus=8
#PBS -l walltime=0:10:00
#PBS -j oe
#PBS -q devel
#PBS -W group_list=s1075

###PBS -S /bin/csh
###PBS -N testlnewc
###PBS -o outfile
###PBS -l select=256:ncpus=8
###PBS -l walltime=1:00:00
###PBS -j oe

setenv DATA_DIRECTORY $PBS_O_WORKDIR/data
setenv RESTART_DIRECTORY $PBS_O_WORKDIR/restart_files
setenv INPUT_DIRECTORY $PBS_O_WORKDIR
setenv SOURCE_DIRECTORY $PBS_O_WORKDIR
mkdir -p $DATA_DIRECTORY
mkdir -p $RESTART_DIRECTORY

cd $PBS_O_WORKDIR
lfs setstripe data          -s  0 -i -1 -c -1
lfs setstripe restart_files -s 0 -i -1 -c  1

cp $INPUT_DIRECTORY/finput.run8b.dat $DATA_DIRECTORY/finput.dat
cp ${HOME}/.cleanup_status_init $DATA_DIRECTORY/.cleanup_status


module purge
module load comp/intel
module load mpi-sgi
##module load comp-pgi/13.7
##module load mpi-sgi/mpt.2.06rp16

set verbose
setenv MPI_TYPE_MAX 65536
setenv MPI_REQUEST_MAX 65536

mpiexec ./3dh > 3dhout
