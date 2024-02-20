#!/bin/tcsh
#PBS -l nodes=1:hima:ppn=64
#PBS -l walltime=0:10:00
#PBS -N collect_results
#PBS -k oe
#PBS -m abe
#PBS -M jasimonpietri01@wm.edu
cd $PBS_O_WORKDIR
conda activate openmm_mpi_cuda
python3 MDSim/scripts/collect_results.py MDSim/logs/automate_parameters38.ini
conda deactivate