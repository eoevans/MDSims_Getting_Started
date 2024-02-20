#!/bin/tcsh
#PBS -l nodes=1:hima:ppn=32
#PBS -l walltime=00:10:00
#PBS -N packmol_vmd_test
#PBS -k oe
#PBS -m abe
#PBS -M jasimonpietri01@wm.edu
cd $PBS_O_WORKDIR
conda activate openmm_mpi_cuda
python3 MDSim/scripts/make_box.py MDSim/data/simulation_results/glycerol/alumina/50A/run3/files.ini
conda deactivate