import configparser
import os
import sys
import traceback
import subprocess

def validate_section(file, section_str):
    valid_section = read_config('MDSim/logs/valid_params.ini', section_str)
    section = read_config(file, section_str)
    section_dict = {}
    for key in valid_section:
        if valid_section[key] == 'None':
            section_dict[key] = section[key]
        elif valid_section[key] == 'int':
            section_dict[key] = section.getint(key)
        elif valid_section[key] == 'float':
            section_dict[key] = section.getfloat(key)
        elif valid_section[key] == 'bool':
            section_dict[key] = section.getboolean(key)
        elif valid_section[key] == 'file':
            if os.path.isfile(section[key]):
                section_dict[key] = section[key]
            else:
                raise FileNotFoundError
        elif valid_section[key] == 'node':
            if section[key] != 'bora' and section[key] != 'hima' and section[key] != 'vortex':
                raise Exception('Error in section '+section+'. '+key+' is not bora or hima.')
    return section_dict

# read in .ini file as a ConfigParser object. returns a ConfigParser object if no section is specified.
# otherwise returns a section of the ConfigParser
def read_config(file, section_str=None):
    if section_str == None:
        config = configparser.ConfigParser(inline_comment_prefixes='#')
        config.read(file)
        return config
    else:
        config = configparser.ConfigParser(inline_comment_prefixes='#')
        config.read(file)
        return config[section_str]

# this is called for every script that is run in MDSim_automate.py (except collect_results.py)
# returns the all files in files.ini and params for section, including the shared params
def get_files_params(file, section_str):
    files = read_config(file, 'files')
    return files, {**validate_section(files['script_params'], section_str), **validate_section(files['script_params'], 'shared')} 

# validates all sections in a .ini. returns True if valid and False if invalid
def validate_file(file):
    params = read_config(file)
    for section in params:
        try:
            validate_section(file, section)
        except Exception:
            print('Error in section '+section)
            traceback.print_exc()
            return False
    return True

def get_surface_dims(file):
    dims = [0,0,0]
    if file.lower() == 'none':
        return dims
    with open(os.path.join(file), 'r') as f:
        lines = f.readlines()
    for i in range(len(lines)):
        args = lines[i].split(' = ')
        if args[0] == ' SET A':
            dims[0] = float(args[1])
            dims[1] = float(lines[i+1].split(' = ')[1])
            dims[2] = float(lines[i+2].split(' = ')[1])
            return dims

# there is a small chance that packmol generates an invalid box (meaning openmm raises an error when trying to use it for a simulation)
# for the run started on 6/27/23, this happened 2/35 times
# for an unknown reason, intra_configparser.py will get stuck when encountering an empty trajectory file
# This function is implemented to get around that, so that if intra_configparser.py finds an empty trajectory,
# it will get an error and continue to the other boxes
def check_error(file):
    if os.stat(file).st_size == 0:
        raise Exception(file+' is empty.')
    
def has_key(section, key_to_find):
    for key in section:
        if key == key_to_find:
            return True
    else:
        return False

def new_trajectory_analysis(files, new_dir):
    # creates a new directory directory "new_dir", inside the directory that holds the trajectory
    # also creates a new files.ini with inter, intra, acf_log, results updated to "new_dir"
    # and creates dir/logs, dir/results, and write files.ini
    # will be used with slicer, and to test and compare various inter/intra parameters
    new_dir = files['dcd'].split('trajectory')[0]+new_dir+'/'
    make_single_dir(new_dir)
    make_single_dir(new_dir+'results/')
    make_single_dir(new_dir+'logs/')
    files['inter_acf'] = new_dir+'inter_acf_MPI.npy'
    files['intra_acf'] = new_dir+'intra_acf_MPI.npy'
    files['acf_log'] = new_dir+'logs/'
    files['results'] = new_dir+'results/results.txt'
    write_config(files, 'files', new_dir+'files.ini')
    return new_dir+'files.ini'


def write_config(section, section_str, file):
    with open(file, 'w') as f:
        f.write('['+section_str+']\n')
        for key in section:
            f.write(key+' = '+section[key]+'\n')
    

def make_single_dir(path):
    if os.path.isdir(path) == False:
        os.mkdir(path)

def write_sh(section, args):
    if section['node_type'] == 'hima':
        if has_key(section, 'gpu') and section.getboolean('gpu'):
            pbsL_string = 'hi07'
        else:
            pbsL_string = '1:hima'
        command = 'python3 '+section['script']+' '
        if has_key(section, 'processors_pn') == False:
            section['processors_pn'] = '32'
    else:
        pbsL_string = section['num_nodes']+':'+section['node_type']
        command = 'mpiexec -np '+section['np']+' -ppn '+str(int(section['np'])//int(section['num_nodes']))+' python3 '+section['script']+' '
    with open (section['submit_file'], 'w') as f:
        f.write('#!/bin/tcsh\n')
        f.write('#PBS -l nodes='+pbsL_string+':ppn='+section['processors_pn']+'\n')
        f.write('#PBS -l walltime='+section['job_time']+'\n')
        f.write('#PBS -N '+section['job_name']+'\n')
        f.write('#PBS -k oe\n#PBS -m abe\n#PBS -M jasimonpietri01@wm.edu\ncd $PBS_O_WORKDIR\n')
        f.write('conda activate '+section['env']+'\n')
        for i in args:
            f.write(command+i+'\n')
        f.write('conda deactivate')

def submit(section, job_nums=[None]):
    if job_nums[0] == None:
        job_nums = []
    print(job_nums)
    if section['run'] == 'yes':
        depend = []
        if len(job_nums) != 0:
            depend = ['-W','depend=afterany']
            for i in job_nums:
                depend[1] = depend[1]+':'+i
        tmp = subprocess.run(['qsub']+depend+[section['submit_file']], capture_output=True, text=True).stdout[:-1]
        print(tmp)
        return tmp
    else:
        return

def get_job_number(job_name):
    qstat = subprocess.run(['qstat','-u','jasimonpietri01'], capture_output=True, text=True).stdout.split()
    for i in range(len(qstat)):
        if qstat[i] == job_name:
            return qstat[i-3]

def submit_fast_job(script, file, arg, job_dependency=[None], env='openmm_mpi_cuda'):
    job_name = script.split('/')[-1].split('.')[0]
    submit_dict = {'node_type': 'hima', 'script': script, 'processors_pn': '1', 'job_time': '00:00:10',
                   'job_name': job_name, 'env': env, 'run': 'yes', 'submit_file': file}
    write_sh(submit_dict, [arg])
    submit(submit_dict, job_dependency)

def get_sim_dir(path):
    path = path.split('/')
    sim_dir = ''
    for i in range(len(path)-1):
        sim_dir = sim_dir + path[i] + '/'
    return sim_dir

if __name__ == '__main__':
    with open('MDSim/a.txt', 'w') as f:
        f.write(get_job_number)