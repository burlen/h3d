#!/bin/bash
#PBS -q normal
#PBS -l nodes=16:ppn=8:native:noflash
#PBS -l walltime=01:00:00
#PBS -N insitu
#PBS -o insitu.out
#PBS -e insitu.err
#PBS -A use300
#PBS -M mahidhar@sdsc.edu
#PBS -m abe
#PBS -V
#PBS -v Catalina_maxhops=6
# Start of user commands - comments start with a hash sign (#)
cd /oasis/scratch/mahidhar/temp_project/UH3D/insitu/run3D-2
cp input/finput.dat.hxv data/finput.dat
cp .cleanup_status_init data/.cleanup_status
mpirun_rsh -hostfile $PBS_NODEFILE -np 128 MV2_CPU_BINDING_LEVEL=socket MV2_CPU_BINDING_POLICY=scatter PYTHONHOME=/oasis/scratch/poleary/temp_project/gue998/Python SOURCE_DIRECTORY=/oasis/scratch/mahidhar/temp_project/UH3D/h3d_bw/code/laura/ INPUT_DIRECTORY=/oasis/scratch/mahidhar/temp_project/UH3D/insitu/run3D-2/input RESTART_DIRECTORY=/oasis/scratch/mahidhar/temp_project/UH3D/insitu/run3D-2/restart DATA_DIRECTORY=/oasis/scratch/mahidhar/temp_project/UH3D/insitu/run3D-2/data MV2_NUM_HCAS=1 MV2_IBA_HCA=mlx4_0 ./3dh.dp.insitu
/bin/rm restart/rest*
/bin/rm data/*gda
/bin/rm *png
cp input/finput.dat.hxv.non data/finput.dat
cp .cleanup_status_init data/.cleanup_status
mpirun_rsh -hostfile $PBS_NODEFILE -np 128 PYTHONHOME=/oasis/scratch/poleary/temp_project/gue998/Python SOURCE_DIRECTORY=/oasis/scratch/mahidhar/temp_project/UH3D/h3d_bw/code/laura/ INPUT_DIRECTORY=/oasis/scratch/mahidhar/temp_project/UH3D/insitu/run3D-2/input RESTART_DIRECTORY=/oasis/scratch/mahidhar/temp_project/UH3D/insitu/run3D-2/restart DATA_DIRECTORY=/oasis/scratch/mahidhar/temp_project/UH3D/insitu/run3D-2/data MV2_NUM_HCAS=1 MV2_IBA_HCA=mlx4_0 ./3dh.dp.den
/bin/rm restart/rest*
/bin/rm data/*gda
