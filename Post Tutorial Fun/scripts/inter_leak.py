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
from CalcT12_anisotropic import values_of_interest
import subprocess

def write_by_rank(rank, acf):
    acf_dir = files['acf_log']+'inter_acfs/'
    make_single_dir(acf_dir)
    np.save(acf_dir+'inter_acf'+str(rank)+'.npy', acf)

def combine_acf(size):
    combined_acf = 0
    for i in range(size):
        acf_file = files['acf_log']+'inter_acfs/inter_acf'+str(i)+'.npy'
        while os.path.exists(acf_file) == False:
            time.sleep(60)
        combined_acf += np.load(acf_file)
        print('rank '+str(i)+'added to total')
    return combined_acf

def get_cutoff_distance(box_dims, cutoff_distance):
    # args are box dimensions in A. Returns half the smallest dimension of the box or 20 A per Singer
    if cutoff_distance == None:
        return 10
    for i in box_dims:
        if i/2 < cutoff_distance:
            cutoff_distance = i/2
    return cutoff_distance

# return every pair that comes within a given distance of the target hydrogen during the simulation
def get_close_hydrogens(all_hydrogens, selH, radius):
    # get all intra hydrogens
    to_subtract = all_hydrogens.select_atoms('not index '+str(selH.index)+' and resid '+str(selH.resid))
    # remove intra hydrogens from the possible pairs
    all_hydrogens = all_hydrogens.difference(to_subtract)
    # atom group of all hydrogens which have entered the radius at any time
    near_hydrogens = all_hydrogens.select_atoms('around '+str(radius)+' id '+str(selH.id))
    for ts in u.trajectory:
        # for every frame, find which hydrogens have are in the radius
        new_hydrogens = all_hydrogens.select_atoms('around '+str(radius)+' id '+str(selH.id))
        # add new hydrogens to old hydrogens
        near_hydrogens = near_hydrogens.union(new_hydrogens)
    return near_hydrogens

# make function to find the 3theta term
def sph_harm(point, points):
    # make sure np.absolute is not used with m=1, m=2
    # we don't care about which direction the vector ij (that connects hydrogen i and hydrogen j) points
    # but the direction needs to be internally consistent
    v = points-point
    v = v - np.round(v/box_dims)*box_dims
    r = np.linalg.norm(v,axis=1)
    return [(1/r**3)*((3*(v[:,2]/r)**2)-1), ((v[:,0]+(v[:,1]*1j))*v[:,2])/(r**2*r**3), (v[:,0]+(v[:,1]*1j))**2/(r**2*r**3)]

def calc_single_acf(vals,lagf,N):
    # 3 x num_frames x num_pairs
    fvi = fft.fft(vals,n=2*N,axis=1)
    single_acf = np.real(fft.ifft(fvi*np.conjugate(fvi),n=2*N,axis=1) [:,:N,:])
    d = N - np.arange(N)
    d = d[:,np.newaxis]
    d = np.tile(d,(single_acf.shape[0],1,single_acf.shape[2]))
    single_acf = single_acf / d
    single_acf = single_acf[:,:int(lagf*N)+1,:]
    print(np.shape(single_acf))
    single_acf = np.sum(single_acf,axis=2)
    print(np.shape(single_acf))
    return single_acf

def calc_inter(u, all_hydrogens, hydrogen_indices, indices, n_ts, logfile):
    # start the looped analysis
    # initialize empty acf_hold for distributed calculation
    acf_local = np.zeros((3, int(params['lag_fraction']*n_ts)+1))
    rejected_hydrogens = []
    for i in indices:
        # choose a spin
        selH = all_hydrogens.select_atoms('index '+str(hydrogen_indices[i]))
        selH = selH[0]
        # how many molecular pairs will each H atom have
        pairs_positions = get_close_hydrogens(all_hydrogens, selH, get_cutoff_distance(box_dims, params['cutoff_distance'])) #intermolecular
        # split up the pairs into smaller sections. This is done to limit memory usage with large simulation cells
        num_pairs_sections = (len(pairs_positions)//max_pairs_len)+1
        indices = np.linspace(0,len(pairs_positions),num_pairs_sections+1).astype(int)

        dt_string = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        with open(logfile, 'a') as f:
            f.write("Hydrogen number "+str(hydrogen_indices[i])+" START at "+dt_string+".\n")
            f.write('number of pairs = '+str(len(pairs_positions))+'\n')

        for i in range(num_pairs_sections):
            section_positions = pairs_positions[indices[i]:indices[i+1]]
            pairs_vals = np.zeros([3,n_ts,len(section_positions)], dtype='complex')
            # find 3theta term for each time step
            for ts in u.trajectory:
                pairs_vals[:,ts.frame,:] = sph_harm(selH.position, section_positions.positions)
            temp = calc_single_acf(pairs_vals,params['lag_fraction'],n_ts)
            if np.all(abs(temp[:,0]/temp[:,1]) < 10):
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
    #Comm.Barrier()
    # for a reason outside my knowledge, this line of code causes the program to hang with large boxes.
    # I might ask the HPC people about this and implement a fix that does not compromise my morals. 
    # But for now, let's have each rank of the process write out a .npy and have rank 0 go back to combine and delete them all.
    #acf_final = Comm.reduce(inter_acf_local, op=MPI.SUM, root=0)
    write_by_rank(rank, inter_acf_local)

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

    with open('MDSim/a.txt', 'a') as f:
        f.write(str(rank)+'\n')

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

    files, params = get_files_params(sys.argv[1], 'inter')

    # for some reason, if one of the trajectories is empty, the batch run to calculate many inter acfs will stop at that trajectory
    # this is added to raise an error if the trajectory is empty, such that the remaining inter acfs can be calculated
    check_error(files['dcd'])

    # I ran out of memory when doing the IFFT(FFT). My solution is to chop up the array that is [time_steps, pairs]
    # into smaller arrays, such that not all pairs are calculated at once. Here is the max length of aforementioned array's pairs dimension
    # this is really only a problem with vortex c18a nodes. The current value is for use with 5 processes per c18a node
    max_pairs_len = 1000
    surface_dims = get_surface_dims(files['surface_details'])
    box_dims = np.array([surface_dims[0]+params['x_len'], surface_dims[1]+params['y_len'], surface_dims[2]+params['z_len']])

    #load PDB [1] and DCD [2] file for reading in
    # u = mda.Universe(sys.argv[1],sys.argv[2])  # always start with a Universe
    if rank == 0:
        u = mda.Universe(files['pdb_sim'],files['dcd'])
    else:
        u = None

    u = Comm.bcast(u, root=0)

    with open('MDSim/a.txt', 'a') as f:
        f.write(str(rank)+'\n')

    n_ts = len(u.trajectory)

    #select H atoms on all molecules
    all_hydrogens = u.select_atoms("name H*")
    acf_final = np.zeros((3, int(params['lag_fraction']*n_ts)+1))

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

    with open('MDSim/a.txt', 'a') as f:
        f.write(str(rank)+'\n')

    # logfile per rank
    fileout = files['acf_log']+'inter_MPI_log_'+str(rank)+'.log'
    with open(fileout, 'w') as f:
            f.write("MPI Size: "+str(size)+"\n")
            f.write("Atoms in entire universe: "+str(len(u.atoms))+"\n")
            f.write("Residues (molecules) entire in universe: "+str(len(u.residues))+"\n")
            f.write('Hydrogens in entire universe: '+str(len(all_hydrogens))+'\n')
            f.write('Hydrogens considered in Acf calculation: '+str(len(hydrogens_to_run))+'\n')
            f.write("Number of time steps: "+str(n_ts)+"\n")

    hydrogen_indices = hydrogens_to_run.indices
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    with open(fileout, 'a') as f:
        f.write('Start pass 1 date and time: '+dt_string+"\n")
    broadcast_hydrogens(hydrogen_indices, size, rank)

    if Comm.Get_rank() == 0:
        acf_final = combine_acf(size)/len(hydrogens_to_run)
        print("inter_acf:"+str(acf_final))
        with open(fileout, 'a') as f:
            f.write("Final ACF: "+str(acf_final)+"\n")
        np.save(files['inter_acf'], acf_final)
        with open(files['results'], 'a') as f:
            f.write('inter\n')
            f.write('script run time = '+time.strftime('%H:%M:%S', time.gmtime(time.time()-start))+'\n')
            for key in params:
                f.write(key+' = '+str(params[key])+'\n')
            voi = values_of_interest(acf_final, numpy_arr=True, return_acf=False)
            f.write('Tau (ps) = ' + str(voi[0]) + '\n')
            f.write('dw (kHz) = '+str(voi[1])+'\n')
            f.write('T1 = ' + str(voi[2]) + '\n')
            f.write('T2 = ' + str(voi[3]) + '\n')
            dt_string = datetime.now().strftime('%d/%m/%Y %H:%M:%S')
            f.write('appended on '+dt_string+'\n')
            f.write('\n')
        # self destruct. When run on very large systems, the code will hang after completion
        # I gave up on fixing the problem properly months ago. This line kills the job so the next one can start
        submit_fast_job('MDSim/scripts/self_destruct.py', 'MDSim/submit/self_destruct.sh', get_job_number('inter_leak'))