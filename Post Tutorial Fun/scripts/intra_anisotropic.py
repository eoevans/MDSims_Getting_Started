import sys
import numpy as np
# import pyfftw.intrafaces.numpy_fft as fft
import numpy.fft as fft
import MDAnalysis as mda
from datetime import datetime
from mpi4py import MPI
import time
import os
import configparser
from params_validate import *
from CalcT12_anisotropic import values_of_interest

def sph_harm(point, points):
    v = points-point
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
    single_acf = single_acf[:,0:int(lagf*N)+1,:]
    single_acf = np.sum(single_acf,axis=2)
    return single_acf

def calc_intra(u, all_hydrogens, hydrogen_indices, indices, n_ts, logfile):
    # start the looped analysis
    # initialize empty acf_hold for distributed calculation
    acf_local = np.zeros((3, int(params['lag_fraction']*n_ts)+1))
    rejected_hydrogens = []
    for i in indices:
        # choose a spin
        selH = all_hydrogens.select_atoms('index '+str(hydrogen_indices[i]))
        selH = selH[0]
        # how many molecular pairs will each H atom have
        pairs_positions = all_hydrogens.select_atoms("not index "+str(hydrogen_indices[i])+" and resid "+str(selH.resid)) #intramolecular
        pairs_vals = np.zeros([3,n_ts,len(pairs_positions)], dtype='complex')

        dt_string = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        with open(logfile, 'a') as f:
            f.write("Hydrogen number "+str(hydrogen_indices[i])+" START at "+dt_string+".\n")

        for ts in u.trajectory:
            pairs_vals[:,ts.frame,:] = sph_harm(selH.position, pairs_positions.positions)
        temp = calc_single_acf(pairs_vals,params['lag_fraction'],n_ts)
        if np.all(abs(temp[:,0]/temp[:,1]) < 10):
            acf_local += temp
        else:
            rejected_hydrogens.append(hydrogen_indices[i])
            with open(logfile, 'a') as f:
                f.write('rejected hydrogen '+str(hydrogen_indices[i])+'\n')
            print('rejected hydrogen '+str(hydrogen_indices[i]))
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

    # run intramolecular autocorrelation calculation
    with open(fileout, 'a') as f:
        f.write("Rank "+str(rank)+" START.\n")
    intra_acf_local, rejected_hydrogens = calc_intra(u, all_hydrogens, hydrogen_indices, np.arange(start_index, end_index), n_ts, fileout)
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    with open(fileout, 'a') as f:
        f.write("Rank "+str(rank)+" DONE.\n")
        f.write("Local ACF: "+str(intra_acf_local)+"\n")
        f.write("Complete date and time: "+dt_string+"\n")
    acf_final = Comm.reduce(intra_acf_local, op=MPI.SUM, root=0)
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
    files, params = get_files_params(sys.argv[1], 'intra')

    # for some reason, if one of the trajectories is empty, the batch run to calculate many intra acfs will stop at that trajectory
    # this is added to raise an error if the trajectory is empty, such that the remaining intra acfs can be calculated
    check_error(files['dcd'])

    alpha = 1 #2.025568e-15 #m^6 s^-2, this is the constant in front of the whole expression for protons
    # alternatively, 2.025568e-15 #A^6 ps^-2
    ts_fraction = 1 #1/fraction of total time vector to sample
    # I ran out of memory when doing the IFFT(FFT). My solution is to chop up the array that is [time_steps, pairs]
    # into smaller arrays, such that not all pairs are calculated at once. Here is the max length of aforementioned array's pairs dimension
    # this is really only a problem with vortex c18a nodes. The current value is for use with 5 processes per c18a node
    max_pairs_len = 1000
    surface_dims = get_surface_dims(files['surface_details'])
    box_dims = np.array([surface_dims[0]+params['x_len'], surface_dims[1]+params['y_len'], surface_dims[2]+params['z_len']])

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

    # logfile per rank
    fileout = files['acf_log']+'intra_MPI_log_'+str(rank)+'.log'
    with open(fileout, 'w') as f:
            f.write("MPI Size: "+str(size)+"\n")
            f.write("Atoms in entire universe: "+str(len(u.atoms))+"\n")
            f.write("Residues (molecules) entire in universe: "+str(len(u.residues))+"\n")
            f.write('Hydrogens in entire universe: '+str(len(all_hydrogens))+'\n')
            f.write('Hydrogens considered in Acf calculation: '+str(len(hydrogens_to_run))+'\n')
            f.write("Number of time steps: "+str(n_ts)+"\n")

    # unresolved error here. Not sure why, but some hydrogens come back with Nans and inf in the array, and I need to resubmit them.
    # this while loop resubmits incorrect hydrogens until we are done or until the tenth resubmition 
    hydrogen_indices = hydrogens_to_run.indices
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
        print("intra_acf:"+str(acf_final))
        with open(fileout, 'a') as f:
            f.write("Final ACF: "+str(acf_final)+"\n")
        np.save(files['intra_acf'], acf_final)
        with open(files['results'], 'a') as f:
            f.write('intra\n')
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
    else:
        pass
