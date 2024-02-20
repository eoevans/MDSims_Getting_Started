from params_validate import *
import os
import sys

def sim_to_box(file):
    index = file.find('sim')
    file = file[:index] + 'box' + file[index+3:]
    return file

# pdb_sim = MDSim/data/box_data/C8H18/CHARMM//run1/C8H18_sim.pdb
def get_run_folder(file):
    arr = file.split('/')
    dir_constant = ''
    for i in range(len(arr)-2):
        dir_constant = dir_constant + arr[i] + '/'
    run_dir = dir_constant + arr[len(arr)-2] + '/'
    box_dir = dir_constant + 'box' + arr[len(arr)-2][3] + '/'
    if os.path.isdir(run_dir) == False:
        os.rename(box_dir, run_dir)
        return
    if len(os.listdir(run_dir)) == 0:
        os.rmdir(run_dir)
        os.rename(box_dir, run_dir)
        return

files = read_config(sys.argv[1], 'files')
get_run_folder(files['pdb_sim'])
if os.path.isdir(files['pdb_sim']) == False:
    os.rename(sim_to_box(files['pdb_sim']),files['pdb_sim'])
    os.rename(sim_to_box(files['psf']),files['psf'])

