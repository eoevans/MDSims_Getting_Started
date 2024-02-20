#!/bin/tcsh
#PBS -l nodes=2:bora:ppn=12
#PBS -l walltime=6:10:00
#PBS -N intra
#PBS -k oe
#PBS -m abe
#PBS -M jasimonpietri01@wm.edu
cd $PBS_O_WORKDIR
conda activate mda
mpiexec -np 40 -ppn 20 python3 MDSim/scripts/intra_anisotropic.py MDSim/data/simulation_results/acetone/control/50A/run1/files.ini
mpiexec -np 40 -ppn 20 python3 MDSim/scripts/intra_anisotropic.py MDSim/data/simulation_results/acetone/control/50A/run2/files.ini
mpiexec -np 40 -ppn 20 python3 MDSim/scripts/intra_anisotropic.py MDSim/data/simulation_results/acetone/control/50A/run3/files.ini
mpiexec -np 40 -ppn 20 python3 MDSim/scripts/intra_anisotropic.py MDSim/data/simulation_results/acetone/control/100A/run1/files.ini
mpiexec -np 40 -ppn 20 python3 MDSim/scripts/intra_anisotropic.py MDSim/data/simulation_results/acetone/control/100A/run2/files.ini
mpiexec -np 40 -ppn 20 python3 MDSim/scripts/intra_anisotropic.py MDSim/data/simulation_results/acetone/control/100A/run3/files.ini
mpiexec -np 40 -ppn 20 python3 MDSim/scripts/intra_anisotropic.py MDSim/data/simulation_results/acetone/control/200A/run1/files.ini
mpiexec -np 40 -ppn 20 python3 MDSim/scripts/intra_anisotropic.py MDSim/data/simulation_results/acetone/control/200A/run2/files.ini
mpiexec -np 40 -ppn 20 python3 MDSim/scripts/intra_anisotropic.py MDSim/data/simulation_results/acetone/control/200A/run3/files.ini
conda deactivate