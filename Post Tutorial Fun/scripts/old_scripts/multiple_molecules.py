import os
import time
import subprocess

def log(job_name, file):
    job = 'start'
    found = True
    while found:
        found = False
        output = subprocess.check_output(['qstat', '-u', 'jasimonpietri01']).split()
        for i in range(len(output)):
            if output[i].decode('ASCII') == job_name:
                found = True
                if job == 'start':
                    job = 'queued'
                    job_number = output[i-3].decode('ASCII')
                    with open('MDSim/logs/auto_log_'+str(log_num)+'.txt', 'a') as f:
                        f.write(job_name+' '+job_number+'\nqueued on '+time.strftime('%c')+'\n')
                if job == 'queued':
                    if output[i+6].decode('ASCII') == 'R':
                        job = 'running'
                        with open('MDSim/logs/auto_log_'+str(log_num)+'.txt', 'a') as f:
                            f.write('running at '+time.strftime('%c')+'\n')
                time.sleep(5)
    with open('MDSim/logs/auto_log_'+str(log_num)+'.txt', 'a') as f:
        f.write('completed on '+time.strftime('%c')+'\n')
    if os.path.isfile(file):
        with open('MDSim/logs/auto_log_'+str(log_num)+'.txt', 'a') as f:
            f.write('Successful execution of above job.\n\n')
        return
    #else:
        #with open('MDSim/logs/auto_log_'+str(log_num)+'.txt', 'a') as f:
        #    f.write('Erroneous execution of above job. Execution terminated.\n\n')
        #raise Exception(job_name+' output file not generated. see '+job_number+' log file for details')
    
# run_time should be in mins
def calc_job_time(run_time):
    total_time = (run_time+3)*mol_num*box_num
    return str(total_time//60)+':'+str(total_time%60)+':00'

# approximate run time for a single run in minutes as an integer
OPENMM_RUN_TIME = 30
INTER_RUN_TIME = 45
INTRA_RUN_TIME = 6
CALCT12_RUN_TIME = 1

dir_before_mol = '/sciclone/home/jasimonpietri01/MDSim/data'
#dir_before_mol = '~/MDSim/data'
mol_num = input('input the number of molecules to run (input \'all\' for every molecule in box_data): ')
if mol_num == 'all':
    mol = os.listdir(dir_before_mol+'/box_data/')
    #for i in range(len(mol_num)):
     #   if os.path.isdir(mol[i]):
      #      pass
       # else:
        #    mol.pop(i)
    mol_num = len(mol)
else:
    mol_num = int(mol_num)
    mol = [''] * mol_num
    for j in range(mol_num):
        mol[j] = input('input the name of molecule '+str(j+1)+' (the name of the folder with the prmtop and inpcrd): ')
        while os.path.isdir(dir_before_mol+'/box_data/'+mol[j]+'/') == False:
            print('Directory not found, try again.\n')
            mol[j] = input('input the molecule name (the name of the folder with the prmtop and inpcrd): ')

boxdir = '/'+input('input the box folder: ')
box_num = int(input('input the number of boxes per molecule: '))
#boo = input('Do you want to use the sim_anneal protocol [y/n]?: ')
#if boo == 'y':
#    sim_type = 'openmm_sim_anneal.py'
#elif boo == 'n':
#    sim_type = 'openmm_cuda.py'
desc = '/'+input('input the name of the directory which will be created to store the trajectory and other results\n'+
             '(should be descriptive/unique) ')

name = [''] * mol_num
for j in range(mol_num):
    for i in os.listdir(dir_before_mol+'/box_data/'+mol[j]+boxdir+'/'):
        if i[-3:] == 'pdb':
            name[j] = i[:3]

# Make any directories that don't exist
for j in mol:
    if os.path.isdir(dir_before_mol+'/simulation_results/'+j+'/') == False:
        os.mkdir(dir_before_mol+'/simulation_results/'+j+'/')
    if os.path.isdir(dir_before_mol+'/simulation_results/'+j+'/'+desc+'/') == False:
        os.mkdir(dir_before_mol+'/simulation_results/'+j+desc+'/')
    for i in range(int(box_num)):
        if os.path.isdir(dir_before_mol+'/simulation_results/'+j+desc+'/run'+str(i+1)+'/') == False:
            os.mkdir(dir_before_mol+'/simulation_results/'+j+desc+'/run'+str(i+1)+'/')
        if os.path.isdir(dir_before_mol+'/simulation_results/'+j+desc+'/run'+str(i+1)+'/results/') == False:
            os.mkdir(dir_before_mol+'/simulation_results/'+j+desc+'/run'+str(i+1)+'/results/')
        if os.path.isdir(dir_before_mol+'/simulation_results/'+j+desc+'/run'+str(i+1)+'/logs/') == False:
            os.mkdir(dir_before_mol+'/simulation_results/'+j+desc+'/run'+str(i+1)+'/logs/')

with open('MDSim/submit/openmm_auto.sh', 'w') as f:
    f.write('#!/bin/tcsh\n#PBS -l nodes=hi07:ppn=64\n#PBS -l walltime='+calc_job_time(OPENMM_RUN_TIME)+'\n#PBS -N openmm_auto\n'+
    '#PBS -k oe\n#PBS -m abe\n#PBS -M jasimonpietri01@wm.edu\ncd $PBS_O_WORKDIR\nconda activate openmm_mpi_cuda')
    for j in range(mol_num):
        for i in range(box_num):
            f.write('\npython3 ~/MDSim/scripts/openmm_sim_anneal.py '+
            '~/MDSim/data/box_data/'+mol[j]+boxdir+'/'+name[j]+'_box_'+str(i+1)+'.prmtop '+
            '~/MDSim/data/box_data/'+mol[j]+boxdir+'/'+name[j]+'_box_'+str(i+1)+'.inpcrd '+
            '~/MDSim/data/simulation_results/'+mol[j]+desc+'/run'+str(i+1)+'/trajectory.dcd '+
            '~/MDSim/data/simulation_results/'+mol[j]+desc+'/run'+str(i+1)+'/'+name[j]+'_box_'+str(i+1)+'.log '+
            '~/MDSim/data/simulation_results/'+mol[j]+desc+'/run'+str(i+1)+'/results/results.txt')
    f.write('\nconda deactivate')

#with open('MDSim/logs/most_recent_log.txt', 'r') as f:
#    log_num = int(f.read()) + 1
#with open('MDSim/logs/most_recent_log.txt', 'w') as f:
#    f.write(str(log_num))

#with open('MDSim/logs/auto_log_'+str(log_num)+'.txt', 'w') as f:
#    f.write('Starting auto run number '+str(log_num)+'. Data are saved in '+desc+' folder.\n\n')

#job_number = subprocess.run(['qsub', 'MDSim/submit/openmm_auto.sh'], capture_output=True, text=True).stdout[:-1]
#log('openmm_auto', dir_before_mol+'/simulation_results/'+mol[mol_num-1]+desc+'/run'+str(box_num)+'/trajectory.dcd')

with open('MDSim/submit/inter_auto.sh', 'w') as f:
    #f.write('#!/bin/tcsh\n#PBS -l nodes=1:bora:ppn=20\n#PBS -l walltime='+calc_job_time(INTER_RUN_TIME)+'\n#PBS -N inter_auto\n'+
    #'#PBS -k oe\n#PBS -m abe\n#PBS -M jasimonpietri01@wm.edu\n#PBS -W depend=afterany:'+job_number+'\ncd $PBS_O_WORKDIR\nconda activate MDenv')
    f.write('#!/bin/tcsh\n#PBS -l nodes=5:bora:ppn=10\n#PBS -l walltime='+calc_job_time(INTER_RUN_TIME)+'\n#PBS -N inter_auto\n'+
    '#PBS -k oe\n#PBS -m abe\n#PBS -M jasimonpietri01@wm.edu\ncd $PBS_O_WORKDIR\nconda activate MDenv')
    for j in range(mol_num):
        for i in range(box_num):
            f.write('\nmpiexec -np 50 -ppn 10 python3 ~/MDSim/scripts/MPI_inter_test.py '+
            '~/MDSim/data/box_data/'+mol[j]+boxdir+'/'+name[j]+'_box_'+str(i+1)+'.pdb '+
            '~/MDSim/data/simulation_results/'+mol[j]+desc+'/run'+str(i+1)+'/trajectory.dcd '+
            '~/MDSim/data/simulation_results/'+mol[j]+desc+'/run'+str(i+1)+'/inter_acf_MPI.npy '+
            '~/MDSim/data/simulation_results/'+mol[j]+desc+'/run'+str(i+1)+'/logs/inter_MPI_log '+
            '~/MDSim/data/simulation_results/'+mol[j]+desc+'/run'+str(i+1)+'/results/results.txt')
    f.write('\nconda deactivate')

with open('MDSim/submit/intra_auto.sh', 'w') as f:
    #f.write('#!/bin/tcsh\n#PBS -l nodes=1:bora:ppn=20\n#PBS -l walltime='+calc_job_time(INTRA_RUN_TIME)+'\n#PBS -N intra_auto\n'+
    #'#PBS -k oe\n#PBS -m abe\n#PBS -M jasimonpietri01@wm.edu\n#PBS -W depend=afterany:'+job_number+'\ncd $PBS_O_WORKDIR\nconda activate MDenv')
    f.write('#!/bin/tcsh\n#PBS -l nodes=1:bora:ppn=20\n#PBS -l walltime='+calc_job_time(INTRA_RUN_TIME)+'\n#PBS -N intra_auto\n'+
    '#PBS -k oe\n#PBS -m abe\n#PBS -M jasimonpietri01@wm.edu\ncd $PBS_O_WORKDIR\nconda activate MDenv')
    for j in range(mol_num):
        for i in range(box_num):
            f.write('\nmpiexec -np 100 -ppn 20 python3 ~/MDSim/scripts/MPI_intra_test.py '+
            '~/MDSim/data/box_data/'+mol[j]+boxdir+'/'+name[j]+'_box_'+str(i+1)+'.pdb '+
            '~/MDSim/data/simulation_results/'+mol[j]+desc+'/run'+str(i+1)+'/trajectory.dcd '+
            '~/MDSim/data/simulation_results/'+mol[j]+desc+'/run'+str(i+1)+'/intra_acf_MPI.npy '+
            '~/MDSim/data/simulation_results/'+mol[j]+desc+'/run'+str(i+1)+'/logs/intra_MPI_log '+
            '~/MDSim/data/simulation_results/'+mol[j]+desc+'/run'+str(i+1)+'/results/results.txt')
    f.write('\nconda deactivate')

#subprocess.run(['qsub', 'MDSim/submit/intra_auto.sh'])
#job_number = subprocess.run(['qsub', 'MDSim/submit/inter_auto.sh'], capture_output=True, text=True).stdout[:-1]

with open('MDSim/submit/calcT12_auto.sh', 'w') as f:
    #f.write('#!/bin/tcsh\n#PBS -l nodes=1:hima:ppn=64\n#PBS -l walltime='+calc_job_time(CALCT12_RUN_TIME)+'\n#PBS -N calcT12_auto\n'+
    #'#PBS -k oe\n#PBS -m abe\n#PBS -M jasimonpietri01@wm.edu\n#PBS -W depend=afterany:'+job_number+'\ncd $PBS_O_WORKDIR\nconda activate openmm_mpi_cuda')
    f.write('#!/bin/tcsh\n#PBS -l nodes=1:hima:ppn=64\n#PBS -l walltime='+calc_job_time(CALCT12_RUN_TIME)+'\n#PBS -N calcT12_auto\n'+
    '#PBS -k oe\n#PBS -m abe\n#PBS -M jasimonpietri01@wm.edu\ncd $PBS_O_WORKDIR\nconda activate openmm_mpi_cuda')
    for j in range(mol_num):
        for i in range(box_num):
            f.write('\npython3 ~/MDSim/scripts/CalcT12.py '+
            '~/MDSim/data/simulation_results/'+mol[j]+desc+'/run'+str(i+1)+'/intra_acf_MPI.npy '+
            '~/MDSim/data/simulation_results/'+mol[j]+desc+'/run'+str(i+1)+'/inter_acf_MPI.npy '+
            '~/MDSim/data/simulation_results/'+mol[j]+desc+'/run'+str(i+1)+'/results/acf_graphs '+
            '~/MDSim/data/simulation_results/'+mol[j]+desc+'/run'+str(i+1)+'/results/results.txt')
    f.write('\nconda deactivate')

#subprocess.run(['qsub', 'MDSim/submit/calcT12_auto.sh'])

#subprocess.call(['qsub', 'MDSim/submit/inter_auto.sh'])
#log('inter_auto', dir_before_mol+'/simulation_results/'+mol[mol_num-1]+desc+'/run'+str(box_num)+'/inter_acf_MPI.npy')

#subprocess.call(['qsub', 'MDSim/submit/intra_auto.sh'])
#log('intra_auto', dir_before_mol+'/simulation_results/'+mol[mol_num-1]+desc+'/run'+str(box_num)+'/intra_acf_MPI.npy')

#subprocess.call(['qsub', 'MDSim/submit/calcT12_auto.sh'])
#log('calcT12_auto.sh', dir_before_mol+'/simulation_results/'+mol[mol_num-1]+desc+'/run'+str(box_num)+'/results/results.txt')

#with open('MDSim/logs/auto_log_'+str(log_num)+'.txt', 'a') as f:
#    f.write('All jobs completed!')