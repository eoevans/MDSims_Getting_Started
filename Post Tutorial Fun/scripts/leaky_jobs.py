# idea. Inter is simply not working with MDSim_automate with many, large boxes. 
# The batch job will just stop at some point, and I don't know why. 
# my temporary plan is to have a list of all files.ini locations to run.
# I want to submit a single job at a time, and hang a 'waiter' job on it. 
# this job will submit the next inter run and hang the next 'waiter' job.

import sys
from params_validate import *
from constants import *
from MDSim_automate import *
import subprocess

def get_wall_time():
    surface_dims = get_surface_dims('MDSim/data/surface_data/aluminum/20x20x10/aluminum.str')
    box_volume = (10**-3)*(script_params.getfloat('shared', 'x_len')*script_params.getfloat('shared', 'y_len')*(script_params.getfloat('shared', 'z_len')+surface_dims[2]))
    sim_time = script_params.getint('openmm', 'sim_steps')*10**-6
    equil_time = script_params.getint('openmm', 'equil_steps')*0.5*10**-6
    walltime = int(get_single_run_time('inter_sub', sim_time, equil_time, box_volume)*1.5) + 10
    return str(walltime//60)+':'+str(walltime%60)+':00'

if __name__ == '__main__':
    inter_jobs_to_run = read_config(sys.argv[1], 'files')
    inter_files = read_config(inter_jobs_to_run[sys.argv[2]], 'files')
    script_params = read_config(inter_files['script_params'])
    inter_wall_time = get_wall_time()
    inter_dict = configparser.ConfigParser()
    inter_dict['inter_leak'] = {'node_type': 'vortex', 'num_nodes': '15', 'processors_pn': '12', 'np': '75', 'env': 'MDenv_vortex',
                                'job_name': 'inter_leak', 'script': 'MDSim/scripts/inter_leak.py', 'run': 'yes',
                                'job_time': inter_wall_time, 'submit_file': 'MDSim/submit/inter_leak.sh'}
    write_sh(inter_dict['inter_leak'], [inter_jobs_to_run[sys.argv[2]]])
    job_num = submit(inter_dict['inter_leak'])
    submit_fast_job('MDSim/scripts/leaky_jobs.py', 'MDSim/submit/waiter_sub.sh', sys.argv[1]+' '+str(int(sys.argv[2])+1), [job_num])
    