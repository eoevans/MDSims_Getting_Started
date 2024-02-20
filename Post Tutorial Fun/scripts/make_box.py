# This script will make take a valid (after adding name and running pdb_tools) 
# and run packmol and vmd autopsf. 

import configparser
import sys
from constants import *
from params_validate import *
import os
import numpy as np

# write the packmol input file. All surfaces are shifted to have only positive coordinates. 
# So the liquids must be shifted out of the way of the surface, which is done below
def write_inp(params,files,z_min,z_max,segid):
    num_mols = count_mols_in_box(files['mol'], [params['x_len'], params['y_len'], z_max-z_min])
    surface_dims = get_surface_dims(files['surface_details'])
    box_start = [0,0,surface_dims[2]+z_min]
    box_end = [params['x_len'], params['y_len'], z_max+surface_dims[2]]

    # write the packmol input file. see https://m3g.github.io/packmol/userguide.shtml for details on the purpose of each line
    with open(files['inp'], 'w') as f:
        f.write('tolerance 2.0\nfiletype pdb\nseed -1\n')
        f.write('output '+files['pdb_box']+'\n')
        f.write('structure '+files['pdb_valid']+'\n')
        f.write('   segid '+segid+'\n')
        f.write('   number '+str(num_mols)+'\n')
        f.write('   resnumbers 0\n')
        f.write('   inside box '+str(box_start[0])+' '+str(box_start[1])+' '+str(box_start[2])+' '+str(box_end[0])+' '+str(box_end[1])+' '+str(box_end[2])+'\n')
        f.write('end structure\n')
    return num_mols

# return the paths to the pdb files and inp files. 
# because of the spacing between columns in pdb files, only 9999 atoms can fit in a single pdb
def get_pdbs(params, files, inp=True):
    MAX_ATOMS_IN_PDB = 9999
    num_mols = count_mols_in_box(files['mol'], [params['x_len'], params['y_len'], params['z_len']])
    num_pdbs = (num_mols//(MAX_ATOMS_IN_PDB+1))+1
    pdb_files = ['']*num_pdbs
    inp_files = ['']*num_pdbs
    if num_pdbs == 1:
        pdb_files[0] = files['pdb_box']
        inp_files[0] = files['inp']
    else:
        for i in range(num_pdbs):
            pdb_files[i] = files['pdb_box'][:-4]+str(i+1)+'.pdb'
            inp_files[i] = files['inp'][:-4]+str(i+1)+'.inp'
    if inp == True:
        return pdb_files, inp_files
    else:
        return pdb_files

def run_packmol(params, files):
    # write the input files and submit them to packmol
    pdb_files, inp_files = get_pdbs(params, files)
    num_pdbs = len(pdb_files)
    z_lens = np.linspace(0, params['z_len'], num_pdbs+1, dtype=int)
    for i in range(num_pdbs):
        new_files = {'inp': inp_files[i], 'pdb_valid': files['pdb_valid'], 'pdb_box': pdb_files[i],
                     'surface_details': files['surface_details'], 'mol': files['mol']}
        write_inp(params, new_files, z_lens[i], z_lens[i+1], letters[i])
        os.system('packmol < '+new_files['inp'])
    return pdb_files

def write_vmd(files, pdb_files):
    # write the liquid only vmd input file. 
    # see https://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-html/node6.html for details on using vmd's psfgen
    with open(files['vmd_box'], 'w') as f:
        f.write('package require psfgen\n')
        f.write('topology '+files['str']+'\n')
        f.write('topology '+files['rtf']+'\n')
        for i in range(len(pdb_files)):
            f.write('segment '+letters[i]+' {pdb '+pdb_files[i]+'}\n')
        for i in range(len(pdb_files)):
            f.write('coordpdb '+pdb_files[i]+' '+letters[i]+'\n')
        f.write('guesscoord\n')
        f.write('writepsf '+files['psf_box'])

def autopsf(files, pdb_files):
    # generate the psf of the liquid only.  
    write_vmd(files, pdb_files)
    os.system('vmd -dispdev text -startup '+files['vmd_box'])

def write_combine_psfs(files, pdb_files):
    with open(files['vmd_combine'], 'w') as f:
        f.write('package require psfgen\n')
        f.write('readpsf '+files['psf_surface']+'\n')
        f.write('readpsf '+files['psf_box']+'\n')
        f.write('coordpdb '+files['pdb_surface']+'\n')
        for i in pdb_files:
            f.write('coordpdb '+i+'\n')
        f.write('writepsf '+files['psf_sim'])

def combine_psfs(files, pdb_files):
    # combine liquid and surface psfs
    write_combine_psfs(files, pdb_files)
    os.system('vmd -dispdev text -startup '+files['vmd_combine'])

if __name__ == '__main__':
    # generate the system
    # this script works with both liquid/surface and liquid-only boxes. 
    # for liquids only, make the pdb and psf necessary for simulation
    # for liquid/surface, make the liquid pdb and mixed psf. the mixed pdb is made by pdb_merge.py
    letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G']
    files, params = get_files_params(sys.argv[1], 'packmol_vmd')
    pdb_files = run_packmol(params, files)
    autopsf(files, pdb_files)
    if files['pdb_surface'].lower != 'none':
        combine_psfs(files, pdb_files)