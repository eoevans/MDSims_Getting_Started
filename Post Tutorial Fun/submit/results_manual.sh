#!/bin/tcsh
#PBS -l nodes=1:vortex:ppn=1
#PBS -l walltime=0:10:00
#PBS -N results_ani
#PBS -k oe
#PBS -m abe
#PBS -M jasimonpietri01@wm.edu
cd $PBS_O_WORKDIR
conda activate MDenv_vortex
python3 MDSim/scripts/collect_results_anisotropic.py MDSim/logs/results_manual/interface_aluminum/
conda deactivate