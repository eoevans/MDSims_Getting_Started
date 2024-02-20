#!/bin/tcsh
#PBS -l nodes=1:hima:ppn=1
#PBS -l walltime=00:00:10
#PBS -N self_destruct
#PBS -k oe
#PBS -m abe
#PBS -M jasimonpietri01@wm.edu
cd $PBS_O_WORKDIR
conda activate openmm_mpi_cuda
python3 MDSim/scripts/self_destruct.py 10989446
conda deactivate