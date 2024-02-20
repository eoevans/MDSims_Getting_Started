#!/bin/tcsh
#PBS -l nodes=1:hima:ppn=32
#PBS -l walltime=00:10:00
#PBS -N rotate
#PBS -k oe
#PBS -m abe
#PBS -M jasimonpietri01@wm.edu
cd $PBS_O_WORKDIR
conda activate mda
python MDSim/scripts/trajectory_rotate.py MDSim/data/simulation_results/C8H18/aluminum/20A/run1/files.ini
python MDSim/scripts/trajectory_rotate.py MDSim/data/simulation_results/C8H18/aluminum/20A/run2/files.ini
python MDSim/scripts/trajectory_rotate.py MDSim/data/simulation_results/C8H18/aluminum/20A/run3/files.ini
conda deactivate