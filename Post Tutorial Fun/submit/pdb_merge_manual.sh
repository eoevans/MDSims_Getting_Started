#!/bin/tcsh
#PBS -l nodes=1:hima:ppn=32
#PBS -l walltime=0:17:00
#PBS -N pdb_merge_interface
#PBS -k oe
#PBS -m abe
#PBS -M jasimonpietri01@wm.edu
cd $PBS_O_WORKDIR
conda activate mda
python3 MDSim/scripts/pdb_merge.py MDSim/data/simulation_results/glycerol/alumina/50A/run3/files.ini
conda deactivate