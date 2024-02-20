import os
import time
import subprocess
import configparser
import shutil
from params_validate import *
import traceback
from constants import *
from copy import deepcopy

# returns the path to the final simulation_results folder based on the automate parameters. 
# a 3D array of shape [num_mols][num_params][num_boxes]
def get_folders(auto, mols, run_folders, folder_type):
    if folder_type == 'box_folder':
        data_folder = 'box_data'
    if folder_type == 'sim_folder':
        data_folder = 'simulation_results'
    files = [None]*len(mols)
    for i in range(len(files)):
        num_runs = [None]*len(run_folders)
        for k in range(len(num_runs)):
            num_boxes = [None]*auto.getint('num_boxes')
            for j in range(len(num_boxes)):
                num_boxes[j] = 'MDSim/data/'+data_folder+'/'+mols[i]+'/'+auto.get(folder_type)+'/'+run_folders[k]+'/run'+str(j+1)+'/'
            num_runs[k] = num_boxes
        files[i] = num_runs
    return files

# returns the mols to run in a list based on automate parameters
# used in collect_results.py.
def get_mols(auto):
    if auto.get('num_mols') == all:
        mol = os.scandir('MDSim/data/box_data/')
    #for i in range(len(params.getint('num_mols'))):
     #   if os.path.isdir(mol[i]) == False:
      #      mol.pop(i)
        auto['num_mols'] = str(len(mol))
    else:
        mol = [''] * auto.getint('num_mols')
        for j in range(auto.getint('num_mols')):
            mol[j] = auto.get('mol'+str(j+1))
    return mol

def get_params_folders(params_folder):
    auto = read_config(params_folder+'auto_params.ini', 'automate')
    data_folders = [[None]*auto.getint('unique_params'),[None]*auto.getint('unique_params')]
    for i in range(auto.getint('unique_params')):
        tmp = read_config(params_folder+'script_params'+str(i+1)+'.ini', 'folders')
        data_folders[0][i] = tmp['box_folder']
        data_folders[1][i] = tmp['sim_folder']
    return data_folders

def make_single_dir(path):
    if os.path.isdir(path) == False:
        os.mkdir(path)

def trunc_path(path, to_remove):
    new_path = ''
    for i in path.split('/')[:to_remove]:
        new_path = new_path + i + '/'
    print(new_path)
    return new_path

# approximates the walltime for each submit file
# we want a walltime that is small in order to avoid long queue times, but not too small that the job doesn't finish
def calc_job_time(section_str, params_folder, automate):
    if automate['surface_dir'].lower() != 'none':
        surface_type = automate['surface_dir'].split('/')[-2]
        surface_type = automate['surface_dir']+surface_type+'.str'
    else:
        surface_type = 'none'
    surface_dims = get_surface_dims(surface_type)
    run_time_sum = 0
    for i in range(automate.getint('unique_params')):
        script_params = read_config(params_folder+'script_params'+str(i+1)+'.ini')
        # box volume including surface z length in nm^3
        box_volume = (10**-3)*(script_params.getfloat('shared', 'x_len')*script_params.getfloat('shared', 'y_len')*(script_params.getfloat('shared', 'z_len')+surface_dims[2]))
        # simulation time in ns
        sim_time = script_params.getint('openmm', 'sim_steps')*10**-6
        # 1/2 equilibration time in ns. I have found that the trend in openmm walltime is more linear when the equilibration time is halved
        equil_time = script_params.getint('openmm', 'equil_steps')*0.5*10**-6
        run_time_sum += get_single_run_time(section_str, sim_time, equil_time, box_volume)*2
    total_time = int((run_time_sum)*automate.getint('num_mols')*automate.getint('num_boxes'))+5
    return str(total_time//60)+':'+str(total_time%60)+':00'

def threeD_to_oneD(sim_folders):
    oneD_array = ['']*len(sim_folders)*len(sim_folders[0])*len(sim_folders[0][0])
    for i in range(len(sim_folders)):
        for j in range(len(sim_folders[0])):
            for k in range(len(sim_folders[0][0])):
                index = i*len(sim_folders[0])*len(sim_folders[0][0])+j*len(sim_folders[0][0])+k
                oneD_array[index] = sim_folders[i][j][k]+'files.ini'
    return oneD_array


def validate_auto(params):
    config = params['automate']
    if config.get('num_mols') == all:
        #mol = os.scandir('MDSim/data/box_data/')
        #for i in range(len(params.getint('num_mols'))):
         #   if os.path.isdir(mol[i]) == False:
          #      mol.pop(i)
        #config['num_mols'] = str(len(mol))
        pass
    else:
        try:
            config.getint('num_mols')
        except ValueError:
            print('num_mols could not be converted to an int. Resubmit automate_parameters'+str(log_num)+'ini')
            return False
        #mol = [''] * config.getint('num_mols')
        #for j in range(config.getint('num_mols')):
        #    mol[j] = config.get('mol'+str(j+1))
        #    if os.path.isdir('MDSim/data/box_data/'+config.get('mol'+str(j+1))+'/') == False:
        #        print('mol'+str(j)+' not found. Resubmit automate_parameters'+str(log_num)+'ini')
        #        return False
    try:
        config.getint('num_boxes')
    except ValueError:
        print('num_boxes could not be converted to an int. Resubmit automate_parameters'+str(log_num)+'ini')         
        return False
    #for j in mol:
    #    if os.path.isdir('MDSim/data/box_data/'+j+'/'+config.get('box_folder')) == False:
    #        print(j+'/box_folder not found. Resubmit automate_parameters'+str(log_num)+'ini')
    #        return False
    for key in config:
        if key[:3] == 'run':
            try:
                config.getboolean(key)
            except ValueError:
                print(key+' is not yes or no. Resubmit automate_parameters'+str(log_num)+'ini')
                return False
    for i in params.sections():
        if i[-3] == 'sub':
            try:
                validate_section(params_folder+auto_params, i)
            except Exception:
                print('Error in section '+i)
                traceback.print_exc()
                return False
    return True

if __name__ == '__main__':
    # the name of the automation parameters file
    auto_params = 'auto_params.ini'
    # the name of the script parameters file
    script_params_file = 'script_params1.ini'

    # check if the user wishes to reuse parameters
    reuse = input('Reuse old parameters? [Y/N] ')
    if reuse.upper() == 'Y':
        # if yes, read in a "runx" folder with both automation and script parameters
        params_folder = input('Enter path to params folder. ')
        if params_folder[-1] != '/':
            params_folder = params_folder + '/'
        params = read_config(params_folder + auto_params)
        automate = params['automate']
    
    elif reuse.upper() == 'N':
        # if no, create new parameters
        # update the parameters file number so that old parameters don't get overwritten
        with open('MDSim/logs/most_recent_log.txt', 'r') as f:
            log_num = int(f.read()) + 1
        with open('MDSim/logs/most_recent_log.txt', 'w') as f:
            f.write(str(log_num))

        # copy last runs parameters into the new folder. This is done so that the comments are still are found in all parameter files
        params_folder_old = 'MDSim/logs/run'+str(log_num-1)+'/'
        params_folder = 'MDSim/logs/run'+str(log_num)+'/'
        make_single_dir(params_folder)
        shutil.copy(params_folder_old+auto_params, params_folder+auto_params)
        shutil.copy(params_folder_old+script_params_file, params_folder+script_params_file)

        # use 'input' to have the program wait until user confirmation that the params are correct
        input('Edit '+params_folder+auto_params+'. Press enter when finished.')
        params = read_config(params_folder+auto_params)
        print(params)
        print(type(params))
        # have the user edit the params as long as the program detects that they are invalid
        while validate_auto(params) == False:
            input('Press enter when finished.')
            params = read_config(params_folder+auto_params)
        automate = params['automate']

        # repeat editing and validation with script parameters
        # something to note, the script parameters are not read in at this step, meaning the user can go back and repair
        # any mistakes before submitting the final script parameter file
        for i in range(automate.getint('unique_params')):
            input('Edit '+'MDSim/logs/run'+str(log_num)+'/script_params_file'+str(i+1)+'.ini. Press enter when finished.')
            while validate_file(params_folder+'script_params'+str(i+1)+'.ini') == False:
                input('Error found, re-edit and press enter when finished.')
            if i+1 != automate.getint('unique_params'):
                shutil.copy(params_folder+'script_params'+str(i+1)+'.ini', params_folder+'script_params'+str(i+2)+'.ini') 

        #shutil.move(params_folder_old, 'MDSim/logs/old_params/')

    # load in all molecules as an array. This is done so that the user has the option to run all molecules, which is mostly a thing of the past now
    mol = get_mols(automate)

    # get the folders where the data for each unique parameter set will be stored. 
    # This is after the automation folder, but before the run1, run2 ... folders. 
    # an array of shape [2, len(unique_params)]
    params_folders = get_params_folders(params_folder)

    # get the paths to all folders that will store the trajectories, files.ini, etc. 
    sim_folders = get_folders(automate, mol, params_folders[1], 'sim_folder')
    # get the paths to all folders that will store the box.pdbs, box.psfs, etc. 
    box_folders = get_folders(automate, mol, params_folders[0], 'box_folder')

    # make folders that don't exist
    for i in range(len(sim_folders)):
        make_single_dir(trunc_path(sim_folders[i][0][0], -4))
        make_single_dir(trunc_path(box_folders[i][0][0], -4))
        make_single_dir(trunc_path(sim_folders[i][0][0], -3))
        make_single_dir(trunc_path(box_folders[i][0][0], -3))
        for j in range(len(sim_folders[i])):
            make_single_dir(trunc_path(sim_folders[i][j][0], -2))
            make_single_dir(trunc_path(box_folders[i][j][0], -2))
            for k in range(len(sim_folders[i][j])):
                make_single_dir(box_folders[i][j][k])
                make_single_dir(sim_folders[i][j][k])
                make_single_dir(sim_folders[i][j][k]+'results/')
                make_single_dir(sim_folders[i][j][k]+'logs/')

# write a file (called files.ini) which contains the path to all necessary files/folders for all relevant programs
    for i in range(len(sim_folders)):
        for j in range(len(sim_folders[i])):
            for k in range(len(sim_folders[i][j])):
                with open(sim_folders[i][j][k]+'files.ini', 'w') as f:
                    f.write('[files]')
                    f.write('\nmol = '+mol[i])
                    f.write('\nscript_params = '+params_folder+'script_params'+str(j+1)+'.ini')
                    f.write('\nauto_params = '+params_folder+auto_params)
                    f.write('\npdb_valid = MDSim/data/single_molecule_data/'+mol[i]+'/'+mol[i]+'_valid.pdb')
                    f.write('\ninp = '+box_folders[i][j][k]+mol[i]+'_packmol.inp')
                    if automate['surface_dir'].lower() == 'none':
                        f.write('\npdb_surface = none')
                        f.write('\npsf_surface = none')
                        f.write('\nsurface_details = none')
                        f.write('\npdb_box = '+box_folders[i][j][k]+mol[i]+'_sim.pdb')
                        f.write('\npsf_box = '+box_folders[i][j][k]+mol[i]+'_sim.psf')
                    else:
                        surface_type = automate['surface_dir'].split('/')[-2]
                        f.write('\npdb_surface = '+automate['surface_dir']+surface_type+'_valid.pdb')
                        f.write('\npsf_surface = '+automate['surface_dir']+surface_type+'_valid.psf')
                        f.write('\nsurface_details = '+automate['surface_dir']+surface_type+'.str')
                        f.write('\npdb_box = '+box_folders[i][j][k]+mol[i]+'_box.pdb')
                        f.write('\npsf_box = '+box_folders[i][j][k]+mol[i]+'_box.psf')
                    f.write('\nvmd_box = '+box_folders[i][j][k]+'psfgen.vmd')
                    f.write('\nvmd_combine = '+box_folders[i][j][k]+'psf_combine.vmd')
                    f.write('\npsf_sim = '+box_folders[i][j][k]+mol[i]+'_sim.psf')
                    f.write('\npdb_sim = '+box_folders[i][j][k]+mol[i]+'_sim.pdb')
                    f.write('\ndcd = '+sim_folders[i][j][k]+'trajectory.dcd')
                    f.write('\ninter_acf = '+sim_folders[i][j][k]+'inter_acf_MPI.npy')
                    f.write('\nintra_acf = '+sim_folders[i][j][k]+'intra_acf_MPI.npy')
                    f.write('\nopenmm_log = '+sim_folders[i][j][k]+'openmm.log')
                    f.write('\nacf_log = '+sim_folders[i][j][k]+'logs/')
                    f.write('\nresults = '+sim_folders[i][j][k]+'results/results.txt')
                    f.write('\ndiffusion_graph = '+sim_folders[i][j][k]+'results/diffusion_graph')
                    f.write(get_CHARMM_files(mol[i]))

    # write every .sh file
    params_files = threeD_to_oneD(sim_folders)
    print(params_files)
    for section in params.sections():
        if section[-3:] == 'sub':
            params[section]['submit_file'] = 'MDSim/submit/'+section+'.sh'
            params[section]['job_time'] = calc_job_time(section, params_folder, automate)
            write_sh(params[section], params_files)

    # submit (or don't) with correct dependancies
    packmol_vmd_job_num = submit(params['packmol_vmd_sub'])
    if automate['surface_dir'].lower() != 'none' and params.getboolean('packmol_vmd_sub', 'run'):
        params['pdb_merge_sub']['run'] = 'yes'
        packmol_vmd_job_num = submit(params['pdb_merge_sub'], [packmol_vmd_job_num])
    openmm_job_num = submit(params['openmm_sub'],[packmol_vmd_job_num])
    inter_job_num = submit(params['inter_sub'],[openmm_job_num])
    intra_job_num = submit(params['intra_sub'],[openmm_job_num])
    diffusion_job_num = submit(params['diffusion_sub'],[openmm_job_num])
    make_single_dir('MDSim/data/numbers_and_graphs/'+automate.get('sim_folder')+'/')
    results_job_num = submit(params['results_sub'],[intra_job_num,inter_job_num])