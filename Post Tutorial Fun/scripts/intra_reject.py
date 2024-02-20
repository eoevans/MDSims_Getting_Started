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

# make function to find the 3theta term
def three_theta(point, points):
    v = points-point
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

def calc_intra(u, all_hydrogens, hydrogen_indices, indeces, n_ts, logfile):
    # start the looped analysis
    # initialize empty acf_hold for distributed calculation
    acf_local = 0
    rejected_hydrogens = []
    for i in indeces:
        # choose a spin
        selH = all_hydrogens.select_atoms('index '+str(hydrogen_indices[i]))
        selH = selH[0]
        # how many molecular pairs will each H atom have
        pairs_positions = all_hydrogens.select_atoms("not index "+str(hydrogen_indices[i])+" and resid "+str(selH.resid)) #intramolecular
        # initialize empty array that is (time steps, connections)
        pairs_vals = np.empty([n_ts,len(pairs_positions)])

        dt_string = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        with open(logfile, 'a') as f:
            f.write("Hydrogen number "+str(hydrogen_indices[i])+" START at "+dt_string+".\n")
        
        # find 3theta term for each time step
        for ts in u.trajectory:
            pairs_vals[ts.frame,:] = three_theta(selH.position, pairs_positions.positions)
        temp = calc_single_acf(pairs_vals,params['lag_fraction'])
        if abs(temp[0]/temp[1]) < 10:
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

    # run intermolecular autocorrelation calculation
    with open(fileout, 'a') as f:
        f.write("Rank "+str(rank)+" START.\n")
    inter_acf_local, rejected_hydrogens = calc_intra(u, all_hydrogens, hydrogen_indices, np.arange(start_index, end_index), n_ts, fileout)
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
    

start = time.time()
files, params = get_files_params(sys.argv[1], 'intra')

# for some reason, if one of the trajectories is empty, the batch run to calculate many intra acfs will stop at that trajectory
# this is added to raise an error if the trajectory is empty, such that the remaining intra acfs can be calculated
check_error(files['dcd'])

alpha = 1 #2.025568e-15 #m^6 s^-2, this is the constant in front of the whole expression for protons
# alternatively, 2.025568e-15 #A^6 ps^-2
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
all_hydrogens = u.select_atoms("name H*")
acf_final = np.zeros(int(params['lag_fraction']*n_ts)+1)
# logfile per rank
fileout = files['acf_log']+'intra_MPI_log_'+str(rank)+'.log'
with open(fileout, 'w') as f:
        f.write("MPI Size: "+str(size)+"\n")
        f.write("Atoms in entire universe: "+str(len(u.atoms))+"\n")
        f.write("Residues (molecules) entire in universe: "+str(len(u.residues))+"\n")
        f.write('Hydrogens in entire universe: '+str(len(all_hydrogens))+'\n')
        f.write('Hydrogens considered in Acf calculation: '+str(len(all_hydrogens))+'\n')
        f.write("Number of time steps: "+str(n_ts)+"\n")

hydrogen_indices = all_hydrogens.indices
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
    acf_final = acf_final/len(all_hydrogens)
    print("intra_acf:"+str(acf_final))
    with open(fileout, 'a') as f:
        f.write("Final ACF: "+str(acf_final)+"\n")
    np.save(files['intra_acf'], acf_final)
    with open(files['results'], 'a') as f:
        f.write('intra\n')
        f.write('script run time = '+time.strftime('%H:%M:%S', time.gmtime(time.time()-start))+'\n')
        for key in params:
            f.write(key+' = '+str(params[key])+'\n')
        f.write('\n')
else:
    pass
