#!/bin/tcsh
#PBS -l nodes=hi07:ppn=32
#PBS -l walltime=1:00:00
#PBS -N openmm
#PBS -k oe
#PBS -m abe
#PBS -M jasimonpietri01@wm.edu
cd $PBS_O_WORKDIR
conda activate openmm_mpi_cuda
python3 MDSim/scripts/openmm_CHARMM.py MDSim/data/simulation_results/water/graphene/20A/run1/files.ini
conda deactivate