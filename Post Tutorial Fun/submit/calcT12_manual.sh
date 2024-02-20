#!/bin/tcsh
#PBS -l nodes=1:hima:ppn=1
#PBS -l walltime=0:10:00
#PBS -N CalcT12_rotated
#PBS -k oe
#PBS -m abe
#PBS -M jasimonpietri01@wm.edu
cd $PBS_O_WORKDIR
conda activate openmm_mpi_cuda
python3 MDSim/scripts/CalcT12_configparser.py MDSim/data/simulation_results/C8H18/CHARMM/run1/cutoff_15A/files.ini
python3 MDSim/scripts/CalcT12_configparser.py MDSim/data/simulation_results/C8H18/CHARMM/run1/cutoff_20A/files.ini
conda deactivate