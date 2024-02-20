import sys
import numpy as np
# import pyfftw.interfaces.numpy_fft as fft
import numpy.fft as fft
import MDAnalysis as mda
from datetime import datetime
from mpi4py import MPI
import time
import os
import configparser
from params_validate import *
from copy import deepcopy
from CalcT12_configparser import values_of_interest

# inter_slice_fast.py functions differently from intra_slice.py. This script identifies the hydrogens inside the slice on the 
# first frame and calculates their acf using FFT. We do not account for whether these hydrogens leave the slice or 
# whether new hydrogens enter. For this, see inter_slice_slow.py

def make_single_dir(path):
    if os.path.isdir(path) == False:
        os.mkdir(path)
    return path

def get_sim_dir(path):
    path = path.split('/')
    sim_dir = ''
    for i in range(len(path)-1):
        sim_dir = sim_dir + path[i] + '/'
    return sim_dir

# make function to find the 3theta term
def three_theta(point, points):
    v = np.absolute(points-point)
    v = v - np.trunc(v/box_dims+0.5)*box_dims
    az = v[:,2]
    a = np.linalg.norm(v,axis=1)
    val_out = (1/a**3)*((3*az**2/a**2)-1)
    return val_out

def calc_single_acf(vals,lagf):
    N = len(vals)
    fvi = fft.rfft(vals,n=2*N,axis=0)
    single_acf = np.real(fft.irfft(fvi*np.conjugate(fvi),n=2*N,axis=0) [:N,:])
    d = N - np.arange(N)
    d = d[:,np.newaxis]
    d = np.tile(d,(1,single_acf.shape[1]))
    single_acf = single_acf / d 
    single_acf = single_acf[0:int(lagf*len(single_acf))+1,:]
    single_acf = np.sum(single_acf,axis=1)
    return single_acf

def calc_inter(u, all_hydrogens, indices, hydrogen_indices, logfile):
    n_ts = len(u.trajectory)
    # start the looped analysis
    # initialize empty acf_hold for distributed calculation
    acf_local = 0
    rejected_hydrogens = []
    for i in indices:
        # choose a spin
        selH = all_hydrogens.select_atoms('index '+str(hydrogen_indices[i]))
        selH = selH[0]
        # how many molecular pairs will each H atom have

        dt_string = datetime.now().strftime('%d/%m/%Y %H:%M:%S')
        with open(logfile, 'a') as f:
            f.write('Hydrogen number '+str(all_hydrogens[i].id)+' START at '+dt_string+'.\n')
        pairs_positions = get_close_hydrogens(all_hydrogens, selH, box_dims[0]/2)
        # initialize empty array that is (time steps, connections)
        pairs_vals = np.zeros([n_ts,len(pairs_positions)])

        dt_string = datetime.now().strftime('%d/%m/%Y %H:%M:%S')
        with open(logfile, 'a') as f:
            f.write('finished getting revelant pairs at '+dt_string+'.\n')
        
        # find 3theta term for each time step
        for ts in u.trajectory:
            pairs_vals[ts.frame,:] = three_theta(selH.position, pairs_positions.positions)
        dt_string = datetime.now().strftime('%d/%m/%Y %H:%M:%S')
        with open(logfile, 'a') as f:
            f.write('finished getting dipolar E at '+dt_string+'.\n')
        acf_local += calc_single_acf(pairs_vals,params['lag_fraction'])
        dt_string = datetime.now().strftime('%d/%m/%Y %H:%M:%S')
        with open(logfile, 'a') as f:
            f.write('finished getting acf at '+dt_string+'.\n')
    
    #acf_hold = acf_local/n_pairs
    return acf_local, rejected_hydrogens

def broadcast_hydrogens(hydrogen_indices, mpi_size, rank):
    min_atoms_per_process = len(hydrogen_indices) // mpi_size
    remainder_atoms = len(hydrogen_indices) % mpi_size
    start_index = rank * min_atoms_per_process
    end_index = (rank + 1) * min_atoms_per_process
    if rank >= mpi_size - remainder_atoms:
        end_index += rank - (mpi_size - remainder_atoms) + 1
        start_index += rank - (mpi_size - remainder_atoms)
    # record current date and time
    now = datetime.now()
    dt_string = now.strftime('%d/%m/%Y %H:%M:%S')

    with open(fileout, 'a') as f:
        f.write('Start date and time: '+dt_string+'\n')
        f.write('Rank of this process: '+str(rank)+'\n')
        f.write('Atoms per process: '+str(end_index-start_index)+'\n')
        #f.write('Atoms to run on this process: '+str(hydrogen_indices[start_index:end_index].ids)+'\n')
        f.write('Start and stop range on this process: '+str(start_index)+' to '+str(end_index)+'\n')

    # run intermolecular autocorrelation calculation
    with open(fileout, 'a') as f:
        f.write('Rank '+str(rank)+' START.\n')
    inter_acf_local, rejected_hydrogens = calc_inter(u, all_hydrogens, np.arange(start_index, end_index), hydrogen_indices, fileout)
    now = datetime.now()
    dt_string = now.strftime('%d/%m/%Y %H:%M:%S')
    with open(fileout, 'a') as f:
        f.write('Rank '+str(rank)+' DONE.\n')
        f.write('Local ACF: '+str(inter_acf_local)+'\n')
        f.write('Complete date and time: '+dt_string+'\n')
    acf_final = Comm.reduce(inter_acf_local, op=MPI.SUM, root=0)
    indices_to_rerun = Comm.gather(rejected_hydrogens, root=0)
    return acf_final, indices_to_rerun

def merge_arr(indices):
    print(indices)
    merged_arr = []
    for i in indices:
        for j in i:
            merged_arr.append(j)
    return merged_arr

# get all the pairs that come within a given distance of the selected hydrogen
# figured it would be easiest to, if the pair is physically relevant at any time, to calculate it for every time
# rather than calculating it for only the times that are relevant. This way, we have continuous data,
# allowing use to use FFT to calculate the acf. 
def get_close_hydrogens(all_hydrogens, selH, radius):
    # get all intra hydrogens
    to_subtract = all_hydrogens.select_atoms('not index '+str(selH.index)+' and resid '+str(selH.resid))
    # remove intra hydrogens from the possible pairs
    all_hydrogens = all_hydrogens.difference(to_subtract)
    # for every frame, find which hydrogens have are in the radius
    new_hydrogens = all_hydrogens.select_atoms('around '+str(radius)+' id '+str(selH.id), updating=True)
    # atom group of all hydrogens which have entered the radius at any time
    near_hydrogens = all_hydrogens.select_atoms('around '+str(radius)+' id '+str(selH.id), updating=False)
    for ts in u.trajectory:
        # combine new hydrogens and old hydrogens without duplicates
        near_hydrogens = near_hydrogens.symmetric_difference(new_hydrogens)
    return near_hydrogens

def get_max_z(u, atom_name, lq_z_len):
    # find the max height of the solid at t=0
    atom_group = u.select_atoms('name '+atom_name)
    max_z = -10
    for i in atom_group:
        if i.position[2] > max_z and i.position[2] < lq_z_len/2:
            max_z = i.position[2]
    if max_z == -10:
        for i in atom_group:
            if i.position[2] > max_z:
                max_z = i.position[2]
    return max_z

def get_included_hydrogens(u, all_hydrogens, z_min, z_max, wrapped):
    if wrapped:
        interface_hydrogens = all_hydrogens.select_atoms('prop z >= '+str(z_min)+' or prop z <= '+str(z_max))
    else:
        interface_hydrogens = all_hydrogens.select_atoms('prop z >= '+str(z_min)+' and prop z <= '+str(z_max))
    return interface_hydrogens.indices

def get_slice_dims(slice_start, slice_length, lq_z_len, box_z_len):
    if slice_start == 'interface':
        z_min = get_max_z(u, 'IAL', lq_z_len)
    elif slice_start == 'bulk':
        top_surface = get_max_z(u, 'IAL', lq_z_len)
        z_min = top_surface + (0.5*lq_z_len) - (0.5*slice_length)
        if z_min > box_z_len:
            z_min = z_min - box_z_len
    else:
        top_surface = get_max_z(u, 'IAL', lq_z_len)
        z_min = top_surface + float(slice_start)
    z_max = z_min + slice_length
    if z_max > box_z_len:
        z_max = z_max - box_z_len
        wrapped = True
    else:
        wrapped = False
    return z_min, z_max, wrapped

def test(z_min, z_max, wrapped):
    if wrapped:
        print('hydrogens above '+str(z_min)+' or below '+str(z_max))
    else:
        print('hydrogens above '+str(z_min)+' and below '+str(z_max))

if __name__ == '__main__':
    start = time.time()
    files, params = get_files_params(sys.argv[1], 'inter')

    # for some reason, if one of the trajectories is empty, the batch run to calculate many intra acfs will stop at that trajectory
    # this is added to raise an error if the trajectory is empty, such that the remaining intra acfs can be calculated
    check_error(files['dcd'])

    Comm = MPI.COMM_WORLD
    rank = Comm.Get_rank()
    size = Comm.Get_size()

    # get hydrogens in the slice. I want to include three cases for the definition of the slice. 
    # case 1. slice_start = interface. The slice is defined as the top of the interface + slice_length
    # case 2. slice_start = bulk. The slice is defined as the middle of the box +/- 0.5*slice_length
    # case 3. slice_start = x. The slice is defined as the slice start (a float input by the user) + slice_length
    # including the interface and bulk case to start will hopefully make things less complicated. 

    params['slice_start'] = sys.argv[2]
    params['slice_length'] = float(sys.argv[3])

    new_dir = get_sim_dir(files['dcd'])+params['slice_start']+'/'
    files['acf_log'] = new_dir+'/logs/'
    files['inter_acf'] = new_dir+'/inter_acf_MPI.npy'
    files['results'] = new_dir+'/results.txt'

    if rank == 0:
        make_single_dir(new_dir)
        make_single_dir(files['acf_log'])

    surface_dims = get_surface_dims(files['surface_details'])
    box_dims = np.array([surface_dims[0]+params['x_len'], surface_dims[1]+params['y_len'], surface_dims[2]+params['z_len']])

    ts_fraction = 1 #1/fraction of total time vector to sample

    #load PDB [1] and DCD [2] file for reading in
    # u = mda.Universe(sys.argv[1],sys.argv[2])  # always start with a Universe
    if rank == 0:
        u = mda.Universe(files['pdb_sim'],files['dcd'])
    else:
        u = None

    u = Comm.bcast(u, root=0)

    n_ts = len(u.trajectory)//ts_fraction

    #select H atoms on all molecules
    all_hydrogens = u.select_atoms('name H*')
    acf_final = 0
    z_min, z_max, wrapped = get_slice_dims(params['slice_start'], params['slice_length'], params['z_len'], box_dims[2])
    test(z_min, z_max, wrapped)
    hydrogen_indices = get_included_hydrogens(u, all_hydrogens, z_min, z_max, wrapped)
    calc_size = len(hydrogen_indices)
# logfile per rank
    fileout = files['acf_log']+'inter_MPI_log_'+str(rank)+'.log'
    with open(fileout, 'w') as f:
        f.write('MPI Size: '+str(size)+'\n')
        f.write('Atoms in entire universe: '+str(len(u.atoms))+'\n')
        f.write('Residues (molecules) entire in universe: '+str(len(u.residues))+'\n')
        f.write('Hydrogens in entire universe: '+str(len(all_hydrogens))+'\n')
        f.write('Hydrogens considered in Acf calculation: '+str(len(hydrogen_indices))+'\n')
        f.write('Number of time steps: '+str(n_ts)+'\n')

    
    count = 0
    while len(hydrogen_indices) != 0:
        count += 1
        now = datetime.now()
        dt_string = now.strftime('%d/%m/%Y %H:%M:%S')
        with open(fileout, 'a') as f:
            f.write('Start pass '+str(count)+' date and time: '+dt_string+'\n')
        temp, hydrogen_indices = broadcast_hydrogens(hydrogen_indices, size, rank)
        if rank == 0:
            acf_final += temp
            hydrogen_indices = merge_arr(hydrogen_indices)
        hydrogen_indices = Comm.bcast(hydrogen_indices, root=0)
        if count == 10:
            break

    if Comm.Get_rank() == 0:
        acf_final = acf_final/calc_size
        print('inter_acf:'+str(acf_final))
        with open(fileout, 'a') as f:
            f.write('Final ACF: '+str(acf_final)+'\n')
        np.save(files['inter_acf'], acf_final)
        tau, dw, t12 = values_of_interest(acf_final)
        with open(files['results'], 'a') as f:
            f.write('inter\n')
            f.write('script run time = '+time.strftime('%H:%M:%S', time.gmtime(time.time()-start))+'\n')
            for key in params:
                f.write(key+' = '+str(params[key])+'\n')
            f.write(f'{tau = }\n')
            f.write(f'{dw = }\n')
            f.write(f'{t12 = }\n')
            f.write('\n')
    else:
        pass
