import MDAnalysis
from params_validate import *
import sys
from params_validate import *
from MDAnalysis import transformations
import numpy as np
import os
import shutil

def make_single_dir(path):
    if os.path.isdir(path) == False:
        os.mkdir(path)

def rotate(ag, ts):
    rotated = transformations.rotate.rotateby(90, direction=[0,1,0], ag=ag)(ts)

def check(atom_group):
    print(np.max(atom_group.positions, axis=0))
    print(np.min(atom_group.positions, axis=0))

def write_config(section, section_str, file):
    with open(file, 'w') as f:
        f.write('['+section_str+']\n')
        for key in section:
            f.write(key+' = '+section[key]+'\n')

if __name__ == '__main__':
    files, params = get_files_params(sys.argv[1], 'inter')
    surface_dims = get_surface_dims(files['surface_details'])
    box_dims = np.array([surface_dims[0]+params['x_len'], surface_dims[1]+params['y_len'], surface_dims[2]+params['z_len']])
    new_dir = files['dcd'].split('trajectory')[0]+'rotated/'
    make_single_dir(new_dir)
    make_single_dir(new_dir+'results/')
    make_single_dir(new_dir+'logs/')

    u = MDAnalysis.Universe(files['pdb_sim'], files['dcd'])

    files['dcd'] = new_dir+'trajectory.dcd'
    files['inter_acf'] = new_dir+'inter_acf_MPI.npy'
    files['intra_acf'] = new_dir+'intra_acf_MPI.npy'
    files['openmm_log'] = new_dir+'openmm.log'
    files['acf_log'] = new_dir+'logs/'
    files['results'] = new_dir+'results/results.txt'

    with MDAnalysis.Writer(files['dcd'], u.atoms.n_atoms) as W:
        for ts in u.trajectory:
            rotate(u.atoms, ts)
            W.write(u.atoms)

    write_config(files, 'files', new_dir+'files.ini')