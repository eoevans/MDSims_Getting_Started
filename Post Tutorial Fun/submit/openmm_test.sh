#!/bin/tcsh
#PBS -l nodes=hi07:ppn=64
#PBS -l walltime=1:20:00
#PBS -N openmm_NVE
#PBS -k oe
#PBS -m abe
#PBS -M jasimonpietri01@wm.edu
cd $PBS_O_WORKDIR
conda activate openmm_mpi_cuda
python3 MDSim/scripts/openmm_combine.py MDSim/data/simulation_results/acetone/surfaces/run1/params.ini MDSim/data/simulation_results/acetone/surfaces/run1/files.ini
conda deactivate