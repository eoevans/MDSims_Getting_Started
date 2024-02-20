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
from intra_slice import *

# make function to find the 3theta term
def three_theta(point, points):
    v = np.absolute(points-point)
    v = v - np.trunc(v/box_dims+0.5)*box_dims
    az = v[:,2]
    a = np.linalg.norm(v,axis=1)
    val_out = (1/a**3)*((3*az**2/a**2)-1)
    return val_out

def calc_single_acf(vals, lagf):
    N = len(vals)
    single_acf = np.zeros(int(lagf*N))
    # special case where lag time = 0
    acf_pairs = vals*vals
    acf_pairs = acf_pairs[np.nonzero(acf_pairs[:,0])[0],:]
    acf_pairs = np.average(acf_pairs, axis=0)
    single_acf[0] = np.sum(acf_pairs)
    for i in range(1,int(lagf*N)):
        # shift and multiply data
        acf_pairs = vals[i:, :]*vals[:-i, :]
        # remove zeros
        acf_pairs = acf_pairs[np.nonzero(acf_pairs[:,0])[0],:]
        if len(acf_pairs) > 0:
            # average and sum
            acf_pairs = np.average(acf_pairs, axis=0)
            single_acf[i] = np.sum(acf_pairs)
    return single_acf

def calc_inter(u, all_hydrogens, indices, hydrogen_indices, entry_frame, logfile):
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
        pairs_positions = get_close_hydrogens(all_hydrogens, selH, box_dims[0])
        # initialize empty array that is (time steps, connections)
        pairs_vals = np.zeros([n_ts,len(pairs_positions)])

        dt_string = datetime.now().strftime('%d/%m/%Y %H:%M:%S')
        with open(logfile, 'a') as f:
            f.write('finished getting revelant pairs at '+dt_string+'.\n')
        
        # find 3theta term for each time step
        for ts in u.trajectory[entry_frame[i]:]:
            if selH.position[2] >= z_mins[ts.frame] and selH.position[2] <= z_maxes[ts.frame]:
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

def broadcast_hydrogens(hydrogen_indices, entry_frame, mpi_size, rank):
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
    inter_acf_local, rejected_hydrogens = calc_inter(u, all_hydrogens, np.arange(start_index, end_index), hydrogen_indices, entry_frame, fileout)
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
# rather than calculating it for only the times that are relevant
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

if __name__ == '__main__':
    start = time.time()
    files, params = get_files_params(sys.argv[1], 'inter')

    # for some reason, if one of the trajectories is empty, the batch run to calculate many intra acfs will stop at that trajectory
    # this is added to raise an error if the trajectory is empty, such that the remaining intra acfs can be calculated
    check_error(files['dcd'])

    new_dir = make_single_dir(get_sim_dir(files['dcd'])+sys.argv[2])
    files['acf_log'] = make_single_dir(new_dir+'/logs/')
    files['intra_acf'] = new_dir+'/intra_acf_MPI.npy'
    files['results'] = new_dir+'/results.txt'

    surface_dims = get_surface_dims(files['surface_details'])
    box_dims = np.array([surface_dims[0]+params['x_len'], surface_dims[1]+params['y_len'], surface_dims[2]+params['z_len']])

    ts_fraction = 1 #1/fraction of total time vector to sample

    Comm = MPI.COMM_WORLD
    rank = Comm.Get_rank()
    size = Comm.Get_size()

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
    # get hydrogens in the slice
    all_solid = u.select_atoms('name IAL')
    z_mins = get_max_z(u, all_solid)
    interface_size = 3.5
    z_maxes = z_mins + interface_size
    hydrogen_indices, entry_frame = get_included_hydrogens(u, all_hydrogens, z_mins, z_maxes)
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
        temp, hydrogen_indices = broadcast_hydrogens(hydrogen_indices, entry_frame, size, rank)
        if rank == 0:
            acf_final += temp
            hydrogen_indices = merge_arr(hydrogen_indices)
        hydrogen_indices = Comm.bcast(hydrogen_indices, root=0)
        if count == 10:
            break

    if Comm.Get_rank() == 0:
        acf_final = acf_final/calc_size
        print('intra_acf:'+str(acf_final))
        with open(fileout, 'a') as f:
            f.write('Final ACF: '+str(acf_final)+'\n')
        np.save(files['intra_acf'], acf_final)
        tau, dw, t12 = values_of_interest(acf_final)
        with open(files['results'], 'a') as f:
            f.write('intra\n')
            f.write('script run time = '+time.strftime('%H:%M:%S', time.gmtime(time.time()-start))+'\n')
            for key in params:
                f.write(key+' = '+str(params[key])+'\n')
            f.write('interface size = '+interface_size+'\n')
            f.write(f'{tau = }\n')
            f.write(f'{dw = }\n')
            f.write(f'{t12 = }\n')
            f.write('\n')
    else:
        pass
