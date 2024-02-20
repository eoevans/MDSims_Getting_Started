#!/bin/tcsh
#PBS -l nodes=1:hima:ppn=1
#PBS -l walltime=00:10:00
#PBS -N test
#PBS -k oe
#PBS -m abe
#PBS -M jasimonpietri01@wm.edu
cd $PBS_O_WORKDIR
python3 MDSim/scripts/test_general.py