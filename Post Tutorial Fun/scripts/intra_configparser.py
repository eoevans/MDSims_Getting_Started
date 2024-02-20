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

def calc_intra(u, all_hydrogens, sub_hydrogens, n_ts, logfile):
    # start the looped analysis
    # initialize empty acf_hold for distributed calculation
    acf_local = 0
    for i in np.arange(len(sub_hydrogens)):
        # choose a spin
        selH = sub_hydrogens[i]
        # how many molecular pairs will each H atom have
        pairs_positions = all_hydrogens.select_atoms("not index "+str(selH.index)+" and resid "+str(selH.resid)) #intramolecular
        # initialize empty array that is (time steps, connections)
        pairs_vals = np.empty([n_ts,len(pairs_positions)])

        dt_string = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        with open(logfile, 'a') as f:
            f.write("Hydrogen number "+str(sub_hydrogens[i].id)+" START at "+dt_string+".\n")
        
        # find 3theta term for each time step
        for ts in u.trajectory[0:n_ts-1]:
            pairs_vals[ts.frame,:] = three_theta(selH.position, pairs_positions.positions)
        temp = calc_single_acf(pairs_vals,params['lag_fraction'])
        if abs(temp[0]/temp[1]) > 10:
            with open(logfile, 'a') as f:
                f.write('rejected hydrogen '+str(selH.id)+'\n')
            print('rejected hydrogen '+str(selH.id))
        else:
            acf_local += temp
    #acf_hold = acf_local/n_pairs
    return acf_local

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
acf_local = np.empty(int(params['lag_fraction']*n_ts)+1)
acf_final = np.empty(int(params['lag_fraction']*n_ts)+1)

# Get the total number of hydrogens and evenly divide it among the processors
atoms_per_process = int(len(all_hydrogens) // size)
hydrogens_to_run = atoms_per_process*size
startIndex = rank * atoms_per_process
endIndex = (rank + 1) * atoms_per_process

# define the atoms that each process will run
sub_hydrogens = all_hydrogens[startIndex:endIndex]

# logfile per rank
fileout = files['acf_log']+'intra_MPI_log_'+str(rank)+'.log'

# record current date and time
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

with open(fileout, 'w') as f:
    f.write("Start date and time: "+dt_string+"\n")
    f.write("MPI Size: "+str(size)+"\n")
    f.write("Rank of this process: "+str(rank)+"\n")
    f.write("Atoms in entire universe: "+str(len(u.atoms))+"\n")
    f.write("Residues (molecules) entire in universe: "+str(len(u.residues))+"\n")
    f.write('Hydrogens in entire universe: '+str(len(all_hydrogens)))
    f.write('Hydrogens considered in Acf calculation: '+str(hydrogens_to_run))
    f.write("Number of time steps: "+str(n_ts)+"\n")
    f.write("Atoms per process: "+str(atoms_per_process)+"\n")
    f.write("Hydrogens to run on this process: "+str(sub_hydrogens.n_atoms)+"\n")
    f.write("Atoms to run on this process: "+str(sub_hydrogens.ids)+"\n")
    f.write("Start and stop range on this process: "+str(startIndex)+" to "+str(endIndex)+"\n")

# run intermolecular autocorrelation calculation
with open(fileout, 'a') as f:
    f.write("Rank "+str(rank)+" START.\n")
inter_acf_local = calc_intra(u, all_hydrogens, sub_hydrogens, n_ts, fileout) # restore this part for other numpy save file arguments , sys.argv[6]+"_"+str(rank)+".npy", sys.argv[7]+"_"+str(rank)+".npy")

now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
with open(fileout, 'a') as f:
    f.write("Rank "+str(rank)+" DONE.\n")
    f.write("Local ACF: "+str(inter_acf_local)+"\n")
    f.write("Complete date and time: "+dt_string+"\n")
# acf_store_MPI = Comm.gather(inter_acf_local,root=0)

acf_final = Comm.reduce(inter_acf_local, op=MPI.SUM, root=0)


if Comm.Get_rank() == 0:
    acf_final = acf_final/hydrogens_to_run
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
