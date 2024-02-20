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

def get_sim_dir(path):
    path = path.split('/')
    sim_dir = ''
    for i in range(len(path)-1):
        sim_dir = sim_dir + path[i] + '/'
    return sim_dir

def make_single_dir(path):
    if os.path.isdir(path) == False:
        os.mkdir(path)
    return path

# make function to find the 3theta term
def three_theta(point, points):
    v = points-point
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

def calc_intra(u, all_hydrogens, indices, hydrogen_indices, entry_frame, logfile, wrapped):
    n_ts = len(u.trajectory)
    # start the looped analysis
    # initialize empty acf_hold for distributed calculation
    acf_local = 0
    rejected_hydrogens = []
    for i in indices:
        # choose a spin
        selH = all_hydrogens.select_atoms('index '+str(hydrogen_indices[i]))
        try:
            selH = selH[0]
        except:
            with open(logfile, 'w') as f:
                f.write(str(hydrogen_indices[i])+'\n')
                f.write(str(len(selH))+'\n')
                atom = u.select_atoms('index '+str(hydrogen_indices[i]))
                f.write(atom)
        # how many molecular pairs will each H atom have
        pairs_positions = all_hydrogens.select_atoms("not index "+str(hydrogen_indices[i])+" and resid "+str(selH.resid)) #intramolecular
        # initialize empty array that is (time steps, connections)
        pairs_vals = np.zeros([n_ts,len(pairs_positions)])

        dt_string = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        with open(logfile, 'a') as f:
            f.write("Hydrogen number "+str(all_hydrogens[i].id)+" START at "+dt_string+".\n")
        
        # find 3theta term for each time step
        for ts in u.trajectory[entry_frame[i]:]:
            if wrapped[ts.frame]:
                if selH.position[2] >= z_mins[ts.frame] or selH.position[2] <= z_maxes[ts.frame]:
                    pairs_vals[ts.frame,:] = three_theta(selH.position, pairs_positions.positions)
            else:
                if selH.position[2] >= z_mins[ts.frame] and selH.position[2] <= z_maxes[ts.frame]:
                    pairs_vals[ts.frame,:] = three_theta(selH.position, pairs_positions.positions)
        temp = calc_single_acf(pairs_vals,params['lag_fraction'])
        #if abs(temp[0]/temp[1]) < 10:
        acf_local += temp
        #else:
        #    rejected_hydrogens.append(hydrogen_indices[i])
        #    with open(logfile, 'a') as f:
        #        f.write('rejected hydrogen '+str(hydrogen_indices[i])+'\n')
        #    print('rejected hydrogen '+str(hydrogen_indices[i]))
    
    #acf_hold = acf_local/n_pairs
    return acf_local, rejected_hydrogens

def broadcast_hydrogens(hydrogen_indices, entry_frame, mpi_size, rank, wrapped):
    min_atoms_per_process = len(hydrogen_indices) // mpi_size
    remainder_atoms = len(hydrogen_indices) % mpi_size
    start_index = rank * min_atoms_per_process
    end_index = (rank + 1) * min_atoms_per_process
    if rank >= mpi_size - remainder_atoms:
        end_index += rank - (mpi_size - remainder_atoms) + 1
        start_index += rank - (mpi_size - remainder_atoms)
    # record current date and time
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

    with open(fileout, 'a') as f:
        f.write("Start date and time: "+dt_string+"\n")
        f.write("Rank of this process: "+str(rank)+"\n")
        f.write("Atoms per process: "+str(end_index-start_index)+"\n")
        #f.write("Atoms to run on this process: "+str(hydrogen_indices[start_index:end_index].ids)+"\n")
        f.write("Start and stop range on this process: "+str(start_index)+" to "+str(end_index)+"\n")

    # run intermolecular autocorrelation calculation
    with open(fileout, 'a') as f:
        f.write("Rank "+str(rank)+" START.\n")
    inter_acf_local, rejected_hydrogens = calc_intra(u, all_hydrogens, np.arange(start_index, end_index), hydrogen_indices, entry_frame, fileout, wrapped)
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    with open(fileout, 'a') as f:
        f.write("Rank "+str(rank)+" DONE.\n")
        f.write("Local ACF: "+str(inter_acf_local)+"\n")
        f.write("Complete date and time: "+dt_string+"\n")
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

def get_included_hydrogens(u, all_hydrogens, z_mins, z_maxes, wrapped):
    found_indices = []
    entry_frame = np.empty(len(all_hydrogens))
    start = 0
    for ts in u.trajectory:
        if wrapped[ts.frame]:
            interface_hydrogens = all_hydrogens.select_atoms('prop z >= '+str(z_mins[ts.frame])+' or prop z <= '+str(z_maxes[ts.frame]))
        else:
            interface_hydrogens = all_hydrogens.select_atoms('prop z >= '+str(z_mins[ts.frame])+' and prop z <= '+str(z_maxes[ts.frame]))
        # identify hydrogens that have yet to enter the box
        new_indices = np.setdiff1d(interface_hydrogens.indices, found_indices)
        # add the new hydrogens to the hydrogens that previously have entered the box
        found_indices = np.concatenate((found_indices, new_indices))
        # store the entry frame of the new hydrogens
        entry_frame[start:start+len(new_indices)] = np.full(len(new_indices), ts.frame)
        start += len(new_indices)
    entry_frame = entry_frame[~np.isnan(entry_frame)]
    return found_indices.astype(int), entry_frame.astype(int)

def get_slice_dims(slice_start, slice_length, lq_z_len, box_z_len):
    if slice_start == 'interface':
        z_min = get_max_z(u, 'IAL', box_z_len)
        z_max = z_min + slice_length
        z_max[z_max > box_z_len] = z_max[z_max > box_z_len] - box_z_len
        wrapped = z_max < z_min
    elif slice_start == 'bulk':
        top_surface = get_max_z(u, 'IAL', box_z_len)
        z_min = top_surface + (0.5*lq_z_len) - (0.5*slice_length)
        z_max = z_min + slice_length
        z_min[z_min > box_z_len] = z_min[z_min > box_z_len] - box_z_len
        z_max[z_max > box_z_len] = z_max[z_max > box_z_len] - box_z_len
        wrapped = np.full((len(z_min)), False)
    else:
        # need to write this option still
        top_surface = get_max_z(u, 'IAL', box_z_len)
        z_min = top_surface + float(slice_start)
    return z_min, z_max, wrapped

def get_max_z(u, atom_name, box_z_len):
    max_z = -10
    atom_group = u.select_atoms('name '+atom_name)
    for i in atom_group:
        # account for case where part of the surface has shifted across the PB
        if i.position[2] > max_z and i.position[2] < box_z_len/2:
            max_z = i.position[2]
            max_atom = i
    # account for case where the entire surface has shifted across the PB
    if max_z == -10:
        for i in atom_group:
            if i.position[2] > max_z:
                max_z = i.position[2]
                max_atom = i
    max_z = np.zeros(len(u.trajectory))
    for ts in u.trajectory:
        max_z[ts.frame] = max_atom.position[2]
    return max_z

def test(z_min, z_max, wrapped):
    for i in range(len(z_min)):
        if wrapped[i]:
            print('hydrogens above '+str(z_min[i])+' or below '+str(z_max[i]))
        else:
            print('hydrogens above '+str(z_min[i])+' and below '+str(z_max[i]))

if __name__ == '__main__':
    start = time.time()
    files, params = get_files_params(sys.argv[1], 'intra')
    
    Comm = MPI.COMM_WORLD
    rank = Comm.Get_rank()
    size = Comm.Get_size()

    # for some reason, if one of the trajectories is empty, the batch run to calculate many intra acfs will stop at that trajectory
    # this is added to raise an error if the trajectory is empty, such that the remaining intra acfs can be calculated
    check_error(files['dcd'])

    surface_dims = get_surface_dims(files['surface_details'])
    box_dims = np.array([surface_dims[0]+params['x_len'], surface_dims[1]+params['y_len'], surface_dims[2]+params['z_len']])

    params['slice_start'] = sys.argv[2]
    params['slice_length'] = float(sys.argv[3])

    #create slice folder inside sim folder
    new_dir = get_sim_dir(files['dcd'])+params['slice_start']+'/'
    files['acf_log'] = new_dir+'/logs/'
    files['intra_acf'] = new_dir+'/intra_acf_MPI.npy'
    files['results'] = new_dir+'/results.txt'

    if rank == 0:
        make_single_dir(new_dir)
        make_single_dir(files['acf_log'])
    
    ts_fraction = 1 #1/fraction of total time vector to sample

    if rank == 0:
        u = mda.Universe(files['pdb_sim'],files['dcd'])
    else:
        u = None

    u = Comm.bcast(u, root=0)

    n_ts = len(u.trajectory)//ts_fraction

    #select H atoms on all molecules
    all_hydrogens = u.select_atoms("name H*")
    acf_final = 0

    z_mins, z_maxes, wrapped = get_slice_dims(params['slice_start'], params['slice_length'], params['z_len'], box_dims[2])
    test(z_mins, z_maxes, wrapped)
    hydrogen_indices, entry_frame = get_included_hydrogens(u, all_hydrogens, z_mins, z_maxes, wrapped)
    calc_size = len(hydrogen_indices)
    # logfile per rank
    fileout = files['acf_log']+'intra_MPI_log_'+str(rank)+'.log'
    with open(fileout, 'w') as f:
        f.write("MPI Size: "+str(size)+"\n")
        f.write("Atoms in entire universe: "+str(len(u.atoms))+"\n")
        f.write("Residues (molecules) entire in universe: "+str(len(u.residues))+"\n")
        f.write('Hydrogens in entire universe: '+str(len(all_hydrogens))+'\n')
        f.write('Hydrogens considered in Acf calculation: '+str(len(hydrogen_indices))+'\n')
        f.write("Number of time steps: "+str(n_ts)+"\n")

    
    count = 0
    while len(hydrogen_indices) != 0:
        count += 1
        now = datetime.now()
        dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
        with open(fileout, 'a') as f:
            f.write('Start pass '+str(count)+' date and time: '+dt_string+"\n")
        temp, hydrogen_indices = broadcast_hydrogens(hydrogen_indices, entry_frame, size, rank, wrapped)
        if rank == 0:
            acf_final += temp
            hydrogen_indices = merge_arr(hydrogen_indices)
        hydrogen_indices = Comm.bcast(hydrogen_indices, root=0)
        if count == 10:
            break

    if Comm.Get_rank() == 0:
        acf_final = acf_final/calc_size
        print("intra_acf:"+str(acf_final))
        with open(fileout, 'a') as f:
            f.write("Final ACF: "+str(acf_final)+"\n")
        np.save(files['intra_acf'], acf_final)
        tau, dw, t12 = values_of_interest(acf_final)
        with open(files['results'], 'a') as f:
            f.write('intra\n')
            f.write('script run time = '+time.strftime('%H:%M:%S', time.gmtime(time.time()-start))+'\n')
            for key in params:
                f.write(key+' = '+str(params[key])+'\n')
            f.write(f'{tau = }\n')
            f.write(f'{dw = }\n')
            f.write(f'{t12 = }\n')
            f.write('\n')
    else:
        pass
