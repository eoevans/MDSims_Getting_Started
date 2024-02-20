#!/bin/tcsh
#PBS -l nodes=1:hima:ppn=64
#PBS -l walltime=0:10:00
#PBS -N CalcT12_CHARMM
#PBS -k oe
#PBS -m abe
#PBS -M jasimonpietri01@wm.edu
cd $PBS_O_WORKDIR
conda activate MDenv
python3 MDSim/scripts/CalcT12_configparser.py MDSim/logs/automate_parameters98.ini MDSim/data/simulation_results/C5H12/CHARMM/run1/files.ini
python3 MDSim/scripts/CalcT12_configparser.py MDSim/logs/automate_parameters98.ini MDSim/data/simulation_results/C5H12/CHARMM/run2/files.ini
python3 MDSim/scripts/CalcT12_configparser.py MDSim/logs/automate_parameters98.ini MDSim/data/simulation_results/C5H12/CHARMM/run3/files.ini
python3 MDSim/scripts/CalcT12_configparser.py MDSim/logs/automate_parameters98.ini MDSim/data/simulation_results/C5H12/CHARMM/run4/files.ini
python3 MDSim/scripts/CalcT12_configparser.py MDSim/logs/automate_parameters98.ini MDSim/data/simulation_results/C5H12/CHARMM/run5/files.ini
python3 MDSim/scripts/CalcT12_configparser.py MDSim/logs/automate_parameters98.ini MDSim/data/simulation_results/C6H14/CHARMM/run1/files.ini
python3 MDSim/scripts/CalcT12_configparser.py MDSim/logs/automate_parameters98.ini MDSim/data/simulation_results/C6H14/CHARMM/run2/files.ini
python3 MDSim/scripts/CalcT12_configparser.py MDSim/logs/automate_parameters98.ini MDSim/data/simulation_results/C6H14/CHARMM/run3/files.ini
python3 MDSim/scripts/CalcT12_configparser.py MDSim/logs/automate_parameters98.ini MDSim/data/simulation_results/C6H14/CHARMM/run4/files.ini
python3 MDSim/scripts/CalcT12_configparser.py MDSim/logs/automate_parameters98.ini MDSim/data/simulation_results/C6H14/CHARMM/run5/files.ini
python3 MDSim/scripts/CalcT12_configparser.py MDSim/logs/automate_parameters98.ini MDSim/data/simulation_results/C7H16/CHARMM/run1/files.ini
python3 MDSim/scripts/CalcT12_configparser.py MDSim/logs/automate_parameters98.ini MDSim/data/simulation_results/C7H16/CHARMM/run2/files.ini
python3 MDSim/scripts/CalcT12_configparser.py MDSim/logs/automate_parameters98.ini MDSim/data/simulation_results/C7H16/CHARMM/run3/files.ini
python3 MDSim/scripts/CalcT12_configparser.py MDSim/logs/automate_parameters98.ini MDSim/data/simulation_results/C7H16/CHARMM/run4/files.ini
python3 MDSim/scripts/CalcT12_configparser.py MDSim/logs/automate_parameters98.ini MDSim/data/simulation_results/C7H16/CHARMM/run5/files.ini
python3 MDSim/scripts/CalcT12_configparser.py MDSim/logs/automate_parameters98.ini MDSim/data/simulation_results/C10H22/CHARMM/run1/files.ini
python3 MDSim/scripts/CalcT12_configparser.py MDSim/logs/automate_parameters98.ini MDSim/data/simulation_results/C10H22/CHARMM/run2/files.ini
python3 MDSim/scripts/CalcT12_configparser.py MDSim/logs/automate_parameters98.ini MDSim/data/simulation_results/C10H22/CHARMM/run3/files.ini
python3 MDSim/scripts/CalcT12_configparser.py MDSim/logs/automate_parameters98.ini MDSim/data/simulation_results/C10H22/CHARMM/run4/files.ini
python3 MDSim/scripts/CalcT12_configparser.py MDSim/logs/automate_parameters98.ini MDSim/data/simulation_results/C10H22/CHARMM/run5/files.ini
python3 MDSim/scripts/CalcT12_configparser.py MDSim/logs/automate_parameters98.ini MDSim/data/simulation_results/C14H30/CHARMM/run1/files.ini
python3 MDSim/scripts/CalcT12_configparser.py MDSim/logs/automate_parameters98.ini MDSim/data/simulation_results/C14H30/CHARMM/run2/files.ini
python3 MDSim/scripts/CalcT12_configparser.py MDSim/logs/automate_parameters98.ini MDSim/data/simulation_results/C14H30/CHARMM/run3/files.ini
python3 MDSim/scripts/CalcT12_configparser.py MDSim/logs/automate_parameters98.ini MDSim/data/simulation_results/C14H30/CHARMM/run4/files.ini
python3 MDSim/scripts/CalcT12_configparser.py MDSim/logs/automate_parameters98.ini MDSim/data/simulation_results/C14H30/CHARMM/run5/files.ini
python3 MDSim/scripts/CalcT12_configparser.py MDSim/logs/automate_parameters98.ini MDSim/data/simulation_results/C16H34/CHARMM/run1/files.ini
python3 MDSim/scripts/CalcT12_configparser.py MDSim/logs/automate_parameters98.ini MDSim/data/simulation_results/C16H34/CHARMM/run2/files.ini
python3 MDSim/scripts/CalcT12_configparser.py MDSim/logs/automate_parameters98.ini MDSim/data/simulation_results/C16H34/CHARMM/run3/files.ini
python3 MDSim/scripts/CalcT12_configparser.py MDSim/logs/automate_parameters98.ini MDSim/data/simulation_results/C16H34/CHARMM/run4/files.ini
python3 MDSim/scripts/CalcT12_configparser.py MDSim/logs/automate_parameters98.ini MDSim/data/simulation_results/C16H34/CHARMM/run5/files.ini
conda deactivate