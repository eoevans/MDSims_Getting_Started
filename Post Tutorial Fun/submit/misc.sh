#!/bin/tcsh
#PBS -l nodes=1:hima:ppn=32
#PBS -l walltime=00:10:00
#PBS -N MISC
#PBS -k oe
#PBS -m abe
#PBS -M jasimonpietri01@wm.edu
cd $PBS_O_WORKDIR
conda activate mda
python3 MDSim/scripts/atom_position_probability.py MDSim/data/simulation_results/water/aluminum/100A/run1/files.ini
python3 MDSim/scripts/atom_position_probability.py MDSim/data/simulation_results/water/graphene/100A/run1/files.ini
conda deactivate