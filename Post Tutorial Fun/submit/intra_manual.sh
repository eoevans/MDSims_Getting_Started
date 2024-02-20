#!/bin/tcsh
#PBS -l nodes=5:bora:ppn=12
#PBS -l walltime=15:00:00
#PBS -N intra_anisotropic
#PBS -k oe
#PBS -m abe
#PBS -M jasimonpietri01@wm.edu
cd $PBS_O_WORKDIR
conda activate mda
mpiexec -np 100 -ppn 10 python3 MDSim/scripts/intra_anisotropic.py MDSim/data/simulation_results/water/alumina/800A/run1/files.ini
mpiexec -np 100 -ppn 10 python3 MDSim/scripts/intra_anisotropic.py MDSim/data/simulation_results/water/alumina/800A/run2/files.ini
mpiexec -np 100 -ppn 10 python3 MDSim/scripts/intra_anisotropic.py MDSim/data/simulation_results/water/alumina/800A/run3/files.ini
conda deactivate