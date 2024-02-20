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
from CalcT12_configparser import write_voi

# make function to find the 3theta term
def three_theta(point, points):
    #v = np.absolute(points-point)
    #v = v - np.trunc(v/box_dims+0.5)*box_dims
    v = points-point
    v = v - np.round(v/box_dims)*box_dims
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

def calc_inter(u, all_hydrogens, hydrogen_indices, indeces, n_ts, logfile):
    # start the looped analysis
    # initialize empty acf_hold for distributed calculation
    acf_local = 0
    rejected_hydrogens = []
    for i in indeces:
        # choose a spin
        selH = all_hydrogens.select_atoms('index '+str(hydrogen_indices[i]))
        selH = selH[0]
        # how many molecular pairs will each H atom have
        pairs_positions = all_hydrogens.select_atoms("not index "+str(hydrogen_indices[i])+" and not resid "+str(selH.resid)) #intermolecular
        # split up the pairs into smaller sections. This is done to limit memory usage with large simulation cells
        num_pairs_sections = (len(pairs_positions)//max_pairs_len)+1
        indices = np.linspace(0,len(pairs_positions),num_pairs_sections+1).astype(int)

        dt_string = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        with open(logfile, 'a') as f:
            f.write("Hydrogen number "+str(hydrogen_indices[i])+" START at "+dt_string+".\n")
        
        for i in range(num_pairs_sections):
            section_positions = pairs_positions[indices[i]:indices[i+1]]
            pairs_vals = np.empty([n_ts,len(section_positions)])
            # find 3theta term for each time step
            for ts in u.trajectory:
                pairs_vals[ts.frame,:] = three_theta(selH.position, section_positions.positions)
            temp = calc_single_acf(pairs_vals,params['lag_fraction'])
            if abs(temp[0]/temp[1]) < 10:
                acf_local += temp
            else:
                rejected_hydrogens.append(hydrogen_indices[i])
                with open(logfile, 'a') as f:
                    f.write('rejected hydrogen '+str(hydrogen_indices[i])+'\n')
                print('rejected hydrogen '+str(hydrogen_indices[i]))
                break
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
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

    with open(fileout, 'a') as f:
        f.write("Start date and time: "+dt_string+"\n")
        f.write("Rank of this process: "+str(rank)+"\n")
        f.write("Atoms per process: "+str(end_index-start_index)+"\n")
        f.write("Atoms to run on this process: "+str(hydrogen_indices[start_index:end_index])+"\n")
        f.write("Start and stop range on this process: "+str(start_index)+" to "+str(end_index)+"\n")

    # run intermolecular autocorrelation calculation
    with open(fileout, 'a') as f:
        f.write("Rank "+str(rank)+" START.\n")
    inter_acf_local, rejected_hydrogens = calc_inter(u, all_hydrogens, hydrogen_indices, np.arange(start_index, end_index), n_ts, fileout)
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
    
if __name__ == '__main__':
    start = time.time()

    Comm = MPI.COMM_WORLD
    rank = Comm.Get_rank()
    size = Comm.Get_size()

    if len(sys.argv) == 3:
        if rank == 0:
            files = read_config(sys.argv[1], 'files')
            path = new_trajectory_analysis(files, sys.argv[2])
            print(path)
        else:
            path = None
        path = Comm.bcast(path, root=0)
    else:
        path = sys.argv[1]

    files, params = get_files_params(path, 'inter')

    # for some reason, if one of the trajectories is empty, the batch run to calculate many inter acfs will stop at that trajectory
    # this is added to raise an error if the trajectory is empty, such that the remaining inter acfs can be calculated
    check_error(files['dcd'])

    ts_fraction = 1 #1/fraction of total time vector to sample
    # I ran out of memory when doing the IFFT(FFT). My solution is to chop up the array that is [time_steps, pairs]
    # into smaller arrays, such that not all pairs are calculated at once. Here is the max length of aforementioned array's pairs dimension
    # this is really only a problem with vortex c18a nodes. The current value is for use with 5 processes per c18a node
    max_pairs_len = 4000
    surface_dims = get_surface_dims(files['surface_details'])
    box_dims = np.array([surface_dims[0]+params['x_len'], surface_dims[1]+params['y_len'], surface_dims[2]+params['z_len']])

    #load PDB [1] and DCD [2] file for reading in
    # u = mda.Universe(sys.argv[1],sys.argv[2])  # always start with a Universe
    if rank == 0:
        u = mda.Universe(files['pdb_sim'],files['dcd'])
    else:
        u = None

    u = Comm.bcast(u, root=0)

    n_ts = len(u.trajectory)//ts_fraction

    #select H atoms on all molecules
    all_hydrogens = u.select_atoms("name H*")
    acf_final = np.zeros(int(params['lag_fraction']*n_ts)+1)

    if rank == 0:
        # randomly sample all hydrogens
        rand_sel = np.random.permutation(len(all_hydrogens))
        hydrogen_indices = np.sort(rand_sel[0:int(len(all_hydrogens)*params['mol_fraction'])])
        print("rand_sel:"+str(rand_sel))
    else:
        hydrogen_indices = None

    # broadcast list of sub-sampled H atoms to whole group
    hydrogen_indices = Comm.bcast(hydrogen_indices, root=0)
    hydrogens_to_run = all_hydrogens[hydrogen_indices]
    hydrogen_indices = hydrogens_to_run.indices

    # logfile per rank
    fileout = files['acf_log']+'inter_MPI_log_'+str(rank)+'.log'
    with open(fileout, 'w') as f:
            f.write("MPI Size: "+str(size)+"\n")
            f.write("Atoms in entire universe: "+str(len(u.atoms))+"\n")
            f.write("Residues (molecules) entire in universe: "+str(len(u.residues))+"\n")
            f.write('Hydrogens in entire universe: '+str(len(all_hydrogens))+'\n')
            f.write('Hydrogens considered in Acf calculation: '+str(len(hydrogens_to_run))+'\n')
            f.write("Number of time steps: "+str(n_ts)+"\n")

    count = 0
    while len(hydrogen_indices) != 0:
        count += 1
        now = datetime.now()
        dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
        with open(fileout, 'a') as f:
            f.write('Start pass '+str(count)+' date and time: '+dt_string+"\n")
        temp, hydrogen_indices = broadcast_hydrogens(hydrogen_indices, size, rank)
        if rank == 0:
            acf_final += temp
            hydrogen_indices = merge_arr(hydrogen_indices)
        hydrogen_indices = Comm.bcast(hydrogen_indices, root=0)
        if count == 10:
            break

    if Comm.Get_rank() == 0:
        acf_final = acf_final/len(hydrogens_to_run)
        print("inter_acf:"+str(acf_final))
        with open(fileout, 'a') as f:
            f.write("Final ACF: "+str(acf_final)+"\n")
        np.save(files['inter_acf'], acf_final)
        with open(files['results'], 'a') as f:
            f.write('inter\n')
            f.write('script run time = '+time.strftime('%H:%M:%S', time.gmtime(time.time()-start))+'\n')
            for key in params:
                f.write(key+' = '+str(params[key])+'\n')
            write_voi(files['results'], acf_final)
            f.write('\n')
    else:
        pass
