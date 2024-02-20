#!/bin/tcsh
#PBS -l nodes=1:hima:ppn=32
#PBS -l walltime=1:00:00
#PBS -N diffusion_CHARMM
#PBS -k oe
#PBS -m abe
#PBS -M jasimonpietri01@wm.edu
cd $PBS_O_WORKDIR
conda activate mda
python3 MDSim/scripts/diffusion.py MDSim/data/simulation_results/water/tip4p2005/run1/files.ini
conda deactivate