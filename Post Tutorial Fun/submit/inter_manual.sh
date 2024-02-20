#!/bin/tcsh
#PBS -l nodes=5:vortex:ppn=12
#PBS -l walltime=03:00:00
#PBS -N inter_anisotropic
#PBS -k oe
#PBS -m abe
#PBS -M jasimonpietri01@wm.edu
cd $PBS_O_WORKDIR
conda activate MDenv_vortex
mpiexec -np 25 -ppn 5 python3 MDSim/scripts/inter_reject.py MDSim/data/simulation_results/water/graphene/100A/run2/files.ini
mpiexec -np 25 -ppn 5 python3 MDSim/scripts/inter_reject.py MDSim/data/simulation_results/water/graphene/100A/run3/files.ini
conda deactivate