#!/bin/tcsh
#PBS -l nodes=1:hima:ppn=1
#PBS -l walltime=00:00:10
#PBS -N leaky_jobs
#PBS -k oe
#PBS -m abe
#PBS -M jasimonpietri01@wm.edu
cd $PBS_O_WORKDIR
conda activate openmm_mpi_cuda
python3 MDSim/scripts/leaky_jobs.py MDSim/logs/inter_files.ini 13
conda deactivate