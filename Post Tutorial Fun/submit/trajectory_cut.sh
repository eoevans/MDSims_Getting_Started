#!/bin/tcsh
#PBS -l nodes=1:hima:ppn=32
#PBS -l walltime=0:12:00
#PBS -N inter_test
#PBS -k oe
#PBS -m abe
#PBS -M jasimonpietri01@wm.edu
cd $PBS_O_WORKDIR
conda activate mda
python3 MDSim/scripts/truncate_trajectory.py MDSim/data/simulation_results/water/alumina/800A/run1/files.ini short 100
conda deactivate