#!/bin/tcsh
#PBS -l nodes=10:vortex:ppn=12
#PBS -l walltime=50:00:00
#PBS -N inter_slice_fast
#PBS -k oe
#PBS -m abe
#PBS -M jasimonpietri01@wm.edu
cd $PBS_O_WORKDIR
conda activate MDenv_vortex
mpiexec -np 50 -ppn 5 python3 MDSim/scripts/inter_slice_fast.py MDSim/data/simulation_results/water/aluminum/50A/run1/files.ini bulk 10
mpiexec -np 50 -ppn 5 python3 MDSim/scripts/inter_slice_fast.py MDSim/data/simulation_results/water/aluminum/50A/run2/files.ini bulk 10
mpiexec -np 50 -ppn 5 python3 MDSim/scripts/inter_slice_fast.py MDSim/data/simulation_results/water/aluminum/50A/run3/files.ini bulk 10
mpiexec -np 50 -ppn 5 python3 MDSim/scripts/inter_slice_fast.py MDSim/data/simulation_results/water/aluminum/100A/run1/files.ini bulk 10
mpiexec -np 50 -ppn 5 python3 MDSim/scripts/inter_slice_fast.py MDSim/data/simulation_results/water/aluminum/100A/run2/files.ini bulk 10
mpiexec -np 50 -ppn 5 python3 MDSim/scripts/inter_slice_fast.py MDSim/data/simulation_results/water/aluminum/100A/run3/files.ini bulk 10
mpiexec -np 50 -ppn 5 python3 MDSim/scripts/inter_slice_fast.py MDSim/data/simulation_results/water/aluminum/200A/run1/files.ini bulk 10
mpiexec -np 50 -ppn 5 python3 MDSim/scripts/inter_slice_fast.py MDSim/data/simulation_results/water/aluminum/200A/run2/files.ini bulk 10
mpiexec -np 50 -ppn 5 python3 MDSim/scripts/inter_slice_fast.py MDSim/data/simulation_results/water/aluminum/200A/run3/files.ini bulk 10
mpiexec -np 50 -ppn 5 python3 MDSim/scripts/inter_slice_fast.py MDSim/data/simulation_results/glycerol/aluminum/50A/run1/files.ini bulk 10
mpiexec -np 50 -ppn 5 python3 MDSim/scripts/inter_slice_fast.py MDSim/data/simulation_results/glycerol/aluminum/50A/run2/files.ini bulk 10
mpiexec -np 50 -ppn 5 python3 MDSim/scripts/inter_slice_fast.py MDSim/data/simulation_results/glycerol/aluminum/50A/run3/files.ini bulk 10
mpiexec -np 50 -ppn 5 python3 MDSim/scripts/inter_slice_fast.py MDSim/data/simulation_results/glycerol/aluminum/100A/run1/files.ini bulk 10
mpiexec -np 50 -ppn 5 python3 MDSim/scripts/inter_slice_fast.py MDSim/data/simulation_results/glycerol/aluminum/100A/run2/files.ini bulk 10
mpiexec -np 50 -ppn 5 python3 MDSim/scripts/inter_slice_fast.py MDSim/data/simulation_results/glycerol/aluminum/100A/run3/files.ini bulk 10
mpiexec -np 50 -ppn 5 python3 MDSim/scripts/inter_slice_fast.py MDSim/data/simulation_results/glycerol/aluminum/200A/run1/files.ini bulk 10
mpiexec -np 50 -ppn 5 python3 MDSim/scripts/inter_slice_fast.py MDSim/data/simulation_results/glycerol/aluminum/200A/run2/files.ini bulk 10
mpiexec -np 50 -ppn 5 python3 MDSim/scripts/inter_slice_fast.py MDSim/data/simulation_results/glycerol/aluminum/200A/run3/files.ini bulk 10
mpiexec -np 50 -ppn 5 python3 MDSim/scripts/inter_slice_fast.py MDSim/data/simulation_results/acetone/aluminum/50A/run1/files.ini bulk 10
mpiexec -np 50 -ppn 5 python3 MDSim/scripts/inter_slice_fast.py MDSim/data/simulation_results/acetone/aluminum/50A/run2/files.ini bulk 10
mpiexec -np 50 -ppn 5 python3 MDSim/scripts/inter_slice_fast.py MDSim/data/simulation_results/acetone/aluminum/50A/run3/files.ini bulk 10
mpiexec -np 50 -ppn 5 python3 MDSim/scripts/inter_slice_fast.py MDSim/data/simulation_results/acetone/aluminum/100A/run1/files.ini bulk 10
mpiexec -np 50 -ppn 5 python3 MDSim/scripts/inter_slice_fast.py MDSim/data/simulation_results/acetone/aluminum/100A/run2/files.ini bulk 10
mpiexec -np 50 -ppn 5 python3 MDSim/scripts/inter_slice_fast.py MDSim/data/simulation_results/acetone/aluminum/100A/run3/files.ini bulk 10
mpiexec -np 50 -ppn 5 python3 MDSim/scripts/inter_slice_fast.py MDSim/data/simulation_results/acetone/aluminum/200A/run1/files.ini bulk 10
mpiexec -np 50 -ppn 5 python3 MDSim/scripts/inter_slice_fast.py MDSim/data/simulation_results/acetone/aluminum/200A/run2/files.ini bulk 10
mpiexec -np 50 -ppn 5 python3 MDSim/scripts/inter_slice_fast.py MDSim/data/simulation_results/acetone/aluminum/200A/run3/files.ini bulk 10
conda deactivate