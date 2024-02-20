#!/bin/tcsh
#PBS -l nodes=1:hima:ppn=32
#PBS -l walltime=0:23:00
#PBS -N results_ani
#PBS -k oe
#PBS -m abe
#PBS -M jasimonpietri01@wm.edu
cd $PBS_O_WORKDIR
conda activate openmm_mpi_cuda
python3 MDSim/scripts/collect_results_anisotropic.py MDSim/data/simulation_results/acetone/control/50A/run1/files.ini
python3 MDSim/scripts/collect_results_anisotropic.py MDSim/data/simulation_results/acetone/control/50A/run2/files.ini
python3 MDSim/scripts/collect_results_anisotropic.py MDSim/data/simulation_results/acetone/control/50A/run3/files.ini
python3 MDSim/scripts/collect_results_anisotropic.py MDSim/data/simulation_results/acetone/control/100A/run1/files.ini
python3 MDSim/scripts/collect_results_anisotropic.py MDSim/data/simulation_results/acetone/control/100A/run2/files.ini
python3 MDSim/scripts/collect_results_anisotropic.py MDSim/data/simulation_results/acetone/control/100A/run3/files.ini
python3 MDSim/scripts/collect_results_anisotropic.py MDSim/data/simulation_results/acetone/control/200A/run1/files.ini
python3 MDSim/scripts/collect_results_anisotropic.py MDSim/data/simulation_results/acetone/control/200A/run2/files.ini
python3 MDSim/scripts/collect_results_anisotropic.py MDSim/data/simulation_results/acetone/control/200A/run3/files.ini
conda deactivate