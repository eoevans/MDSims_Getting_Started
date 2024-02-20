#!/bin/tcsh
#PBS -l nodes=15:vortex:ppn=12
#PBS -l walltime=3:3:00
#PBS -N inter_leak
#PBS -k oe
#PBS -m abe
#PBS -M jasimonpietri01@wm.edu
cd $PBS_O_WORKDIR
conda activate MDenv_vortex
mpiexec -np 75 -ppn 5 python3 MDSim/scripts/inter_leak.py MDSim/data/simulation_results/acetone/control/200A/run3/files.ini
conda deactivate