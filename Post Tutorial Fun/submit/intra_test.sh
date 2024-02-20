#!/bin/tcsh
#PBS -l nodes=7:vortex:ppn=12
#PBS -l walltime=5:00:00
#PBS -N intra_test
#PBS -k oe
#PBS -m abe
#PBS -M jasimonpietri01@wm.edu
cd $PBS_O_WORKDIR
conda activate MDenv_vortex
mpiexec -np 35 -ppn 5 python3 MDSim/scripts/intra_edit.py MDSim/data/simulation_results/water/aluminum/200A/run1/files.ini
conda deactivate