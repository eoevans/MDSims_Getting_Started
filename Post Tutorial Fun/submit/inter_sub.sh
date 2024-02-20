#!/bin/tcsh
#PBS -l nodes=5:vortex:ppn=12
#PBS -l walltime=13:36:00
#PBS -N inter
#PBS -k oe
#PBS -m abe
#PBS -M jasimonpietri01@wm.edu
cd $PBS_O_WORKDIR
conda activate MDenv_vortex
mpiexec -np 25 -ppn 5 python3 MDSim/scripts/inter_anisotropic.py MDSim/data/simulation_results/acetone/control/50A/run1/files.ini
mpiexec -np 25 -ppn 5 python3 MDSim/scripts/inter_anisotropic.py MDSim/data/simulation_results/acetone/control/50A/run2/files.ini
mpiexec -np 25 -ppn 5 python3 MDSim/scripts/inter_anisotropic.py MDSim/data/simulation_results/acetone/control/50A/run3/files.ini
mpiexec -np 25 -ppn 5 python3 MDSim/scripts/inter_anisotropic.py MDSim/data/simulation_results/acetone/control/100A/run1/files.ini
mpiexec -np 25 -ppn 5 python3 MDSim/scripts/inter_anisotropic.py MDSim/data/simulation_results/acetone/control/100A/run2/files.ini
mpiexec -np 25 -ppn 5 python3 MDSim/scripts/inter_anisotropic.py MDSim/data/simulation_results/acetone/control/100A/run3/files.ini
mpiexec -np 25 -ppn 5 python3 MDSim/scripts/inter_anisotropic.py MDSim/data/simulation_results/acetone/control/200A/run1/files.ini
mpiexec -np 25 -ppn 5 python3 MDSim/scripts/inter_anisotropic.py MDSim/data/simulation_results/acetone/control/200A/run2/files.ini
mpiexec -np 25 -ppn 5 python3 MDSim/scripts/inter_anisotropic.py MDSim/data/simulation_results/acetone/control/200A/run3/files.ini
conda deactivate