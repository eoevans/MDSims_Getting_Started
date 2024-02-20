import os
import time
import subprocess

mol = input('input the name of molecule (the name of the folder with the prmtop and inpcrd): ')
while os.path.isdir('/sciclone/home/jasimonpietri01/MDSim/data/box_data/'+mol+'/') == False:
    print('Directory not found, try again')
    mol = input('input the molecule name (the name of the folder with the prmtop and inpcrd): ')
box_name = input('input the full name of the box: ')
boo = input('Do you want to use the sim_anneal protocol [y/n]?: ')
if boo == 'y':
    sim_type = 'openmm_sim_anneal.py'

def check_for_completion(job_name):
    job = 'start'
    while job != 'completed':
        job = 'completed'
        output = subprocess.check_output(['qstat', '-u', 'jasimonpietri01']).split()
        for i in output:
            if i.decode('ASCII') == job_name:
                job = 'in progress'
                time.sleep(5)
                break

if os.path.isdir('/sciclone/home/jasimonpietri01/MDSim/data/simulation_results/'+mol+'/') == False:
    os.mkdir('/sciclone/home/jasimonpietri01/MDSim/data/simulation_results/'+mol+'/')

if os.path.isdir('/sciclone/home/jasimonpietri01/MDSim/data/simulation_results/'+mol+'/sim_anneal/') == False:
    os.mkdir('/sciclone/home/jasimonpietri01/MDSim/data/simulation_results/'+mol+'/sim_anneal/')

for i in range(int(box_num)):
    if os.path.isdir('/sciclone/home/jasimonpietri01/MDSim/data/simulation_results/'+mol+'/sim_anneal/run'+str(i+1)+'/') == False:
        os.mkdir('/sciclone/home/jasimonpietri01/MDSim/data/simulation_results/'+mol+'/sim_anneal/run'+str(i+1)+'/')
    if os.path.isdir('/sciclone/home/jasimonpietri01/MDSim/data/simulation_results/'+mol+'/sim_anneal/run'+str(i+1)+'/results/') == False:
        os.mkdir('/sciclone/home/jasimonpietri01/MDSim/data/simulation_results/'+mol+'/sim_anneal/run'+str(i+1)+'/results/')

with open('MDSim/submit/openmm_auto.sh', 'w') as f:
    f.write('#!/bin/tcsh\n#PBS -l nodes=hi07:ppn=64\n#PBS -l walltime=5:00:00\n#PBS -N openmm_auto\n'+
    '#PBS -k oe\n#PBS -m abe\n#PBS -M jasimonpietri01@wm.edu\ncd $PBS_O_WORKDIR\nconda activate openmm_mpi_cuda')
    for i in range(int(box_num)):
        f.write('\npython3 ~/MDSim/scripts/openmm_sim_anneal.py '+
        '~/MDSim/data/box_data/'+mol+'/'+name+'_box_'+str(i+1)+'.prmtop '+
        '~/MDSim/data/box_data/'+mol+'/'+name+'_box_'+str(i+1)+'.inpcrd '+
        '~/MDSim/data/simulation_results/'+mol+'/sim_anneal/run'+str(i+1)+'/trajectory.dcd '+
        '~/MDSim/data/simulation_results/'+mol+'/sim_anneal/run'+str(i+1)+'/'+name+'_box_'+str(i+1)+'.log '+
        '~/MDSim/data/simulation_results/'+mol+'/sim_anneal/run'+str(i+1)+'/results/results.txt')
    f.write('\nconda deactivate')

for i in range(int(box_num)):
    if os.path.isdir('/sciclone/home/jasimonpietri01/MDSim/data/simulation_results/'+mol+'/sim_anneal/run'+str(i+1)+'/logs/') == False:
        os.mkdir('/sciclone/home/jasimonpietri01/MDSim/data/simulation_results/'+mol+'/sim_anneal/run'+str(i+1)+'/logs/')

with open('MDSim/submit/inter_auto.sh', 'w') as f:
    f.write('#!/bin/tcsh\n#PBS -l nodes=10:bora:ppn=20\n#PBS -l walltime=8:00:00\n#PBS -N inter_auto\n'+
    '#PBS -k oe\n#PBS -m abe\n#PBS -M jasimonpietri01@wm.edu\ncd $PBS_O_WORKDIR\nconda activate MDenv')
    for i in range(int(box_num)):
        f.write('\nmpiexec -np 40 -ppn 4 python3 ~/MDSim/scripts/MPI_inter_test.py '+
        '~/MDSim/data/box_data/'+mol+'/'+name+'_box_'+str(i+1)+'.pdb '+
        '~/MDSim/data/simulation_results/'+mol+'/sim_anneal/run'+str(i+1)+'/trajectory.dcd '+
        '~/MDSim/data/simulation_results/'+mol+'/sim_anneal/run'+str(i+1)+'/inter_acf_MPI.npy '+
        '~/MDSim/data/simulation_results/'+mol+'/sim_anneal/run'+str(i+1)+'/logs/inter_MPI_log '+
        '~/MDSim/data/simulation_results/'+mol+'/sim_anneal/run'+str(i+1)+'/results/results.txt')
    f.write('\nconda deactivate')

with open('MDSim/submit/intra_auto.sh', 'w') as f:
    f.write('#!/bin/tcsh\n#PBS -l nodes=10:bora:ppn=20\n#PBS -l walltime=2:00:00\n#PBS -N intra_auto\n'+
    '#PBS -k oe\n#PBS -m abe\n#PBS -M jasimonpietri01@wm.edu\ncd $PBS_O_WORKDIR\nconda activate MDenv')
    for i in range(int(box_num)):
        f.write('\nmpiexec -np 40 -ppn 4 python3 ~/MDSim/scripts/MPI_intra_test.py '+
        '~/MDSim/data/box_data/'+mol+'/'+name+'_box_'+str(i+1)+'.pdb '+
        '~/MDSim/data/simulation_results/'+mol+'/sim_anneal/run'+str(i+1)+'/trajectory.dcd '+
        '~/MDSim/data/simulation_results/'+mol+'/sim_anneal/run'+str(i+1)+'/intra_acf_MPI.npy '+
        '~/MDSim/data/simulation_results/'+mol+'/sim_anneal/run'+str(i+1)+'/logs/intra_MPI_log '+
        '~/MDSim/data/simulation_results/'+mol+'/sim_anneal/run'+str(i+1)+'/results/results.txt')
    f.write('\nconda deactivate')

with open('MDSim/submit/calcT12_auto.sh', 'w') as f:
    f.write('#!/bin/tcsh\n#PBS -l nodes=1:hima:ppn=64\n#PBS -l walltime=1:00:00\n#PBS -N calcT12_auto\n'+
    '#PBS -k oe\n#PBS -m abe\n#PBS -M jasimonpietri01@wm.edu\ncd $PBS_O_WORKDIR\nconda activate openmm_mpi_cuda')
    for i in range(int(box_num)):
        f.write('\npython3 ~/MDSim/scripts/CalcT12.py '+
        '~/MDSim/data/simulation_results/'+mol+'/sim_anneal/run'+str(i+1)+'/intra_acf_MPI.npy '+
        '~/MDSim/data/simulation_results/'+mol+'/sim_anneal/run'+str(i+1)+'/inter_acf_MPI.npy '+
        '~/MDSim/data/simulation_results/'+mol+'/sim_anneal/run'+str(i+1)+'/results/acf_graphs '+
        '~/MDSim/data/simulation_results/'+mol+'/sim_anneal/run'+str(i+1)+'/results/results.txt')
    f.write('\nconda deactivate')

subprocess.call(['qsub', 'MDSim/submit/openmm_auto.sh'])
check_for_completion('openmm_auto')

subprocess.call(['qsub', 'MDSim/submit/inter_auto.sh'])
check_for_completion('inter_auto')

subprocess.call(['qsub', 'MDSim/submit/intra_auto.sh'])
check_for_completion('intra_auto')

subprocess.call(['qsub', 'MDSim/submit/calcT12_auto.sh'])
check_for_completion('calcT12_auto.sh')
print('Done!')