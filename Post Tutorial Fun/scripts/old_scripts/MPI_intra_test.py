import sys
import numpy as np
# import pyfftw.interfaces.numpy_fft as fft
import numpy.fft as fft
import MDAnalysis as mda
from datetime import datetime
from mpi4py import MPI
import time
import os

start = time.time()

# make function to find the 3theta term
def three_theta(point, points):
    v = points-point
    az = v[:,2]
    a = np.linalg.norm(v,axis=1)
    val_out = (1/a**3)*((3*az**2/a**2)-1)
    return val_out

def calc_multi_acf(vals,lagf):
    multi_acf = 0
    for i in range(vals.shape[1]):
        N = len(vals)
        fvi = fft.rfft(vals[:,i,:],n=2*N,axis=0)
        single_acf = np.real(fft.irfft(fvi*np.conjugate(fvi),n=2*N,axis=0) [:N,:])
        d = N - np.arange(N)
        d = d[:,np.newaxis]
        d = np.tile(d,(1,single_acf.shape[1]))
        single_acf = single_acf / d 
        single_acf = single_acf[0:int(lagf*len(single_acf))+1,:]
        single_acf = np.sum(single_acf,axis=1)
        multi_acf += single_acf
    return multi_acf

def calc_intra(u, all_hydrogens, sub_hydrogens, n_ts, lag_fraction, logfile):
    # how many molecular pairs will each H atom have
    n_pairs = len(all_hydrogens.select_atoms("resid 1"))-1 #for intramolecular pairs

    # initialize empty arrays that are (time steps, selected atoms, connections)
    # we can do this in one array for intramolecular pairs because there are relatively few of them per molecule
    pairs_vals = np.empty([n_ts,sub_hydrogens.n_atoms,n_pairs])

    # start the looped analysis
    # initialize empty acf_hold for distributed calculation
    acf_hold = 0

    # start the looped analysis
    # choose a spin
    for i in np.arange(len(sub_hydrogens)):
        selH = sub_hydrogens[i]
        dt_string = datetime.now().strftime("%d/%m/%Y %H:%M:%S")

        with open(logfile, 'a') as f:
            f.write("Hydrogen number "+str(sub_hydrogens[i].id)+" START at "+dt_string+".\n")

        # subset other H atoms
        pairs_positions = all_hydrogens.select_atoms("not index "+str(selH.index)+" and resid "+str(selH.resid)) #intramolecular

        # find 3theta term for each time step
        for ts in u.trajectory[0:n_ts-1]:
            pairs_vals[ts.frame,i,:] = three_theta(selH.position, pairs_positions.positions)
    with open(logfile, 'a') as f:
        f.write("pairs_vals.shape: "+str(pairs_vals.shape)+"\n")
        f.write("pairs_vals: "+str(pairs_vals)+"\n")
    acf_hold = calc_multi_acf(pairs_vals,lag_fraction)
    #acf_hold = acf_hold/acf_hold[0]
    return acf_hold

alpha = 1 #2.025568e-15 #m^6 s^-2, this is the constant in front of the whole expression for protons
# alternatively, 2.025568e-15 #A^6 ps^-2
lag_fraction = 0.5  #fraction of whole time vector for lag times
ts_fraction = 1 #1/fraction of total time vector to sample
mol_fraction = 1 #because we sample the entire intramolecular space

Comm = MPI.COMM_WORLD
rank = Comm.Get_rank()
size = Comm.Get_size()

#load PDB [1] and DCD [2] file for reading in
# u = mda.Universe(sys.argv[1],sys.argv[2])  # always start with a Universe
if rank == 0:
    u = mda.Universe(sys.argv[1],sys.argv[2])
else:
    u = None

u = Comm.bcast(u, root=0)

n_ts = len(u.trajectory)//ts_fraction

#select H atoms on all molecules
all_hydrogens = u.select_atoms("name H*")
acf_local = np.empty(int(lag_fraction*n_ts)+1)
acf_final = np.empty(int(lag_fraction*n_ts)+1)

# Get the total number of hydrogens and divide it among the processors
atoms_per_process = len(all_hydrogens) // size
rank = Comm.Get_rank()
startIndex = rank * atoms_per_process
endIndex = (rank + 1) * atoms_per_process

# make a selection string for each process to find its own H atoms
sel_str = '('
for i in range(startIndex,endIndex):
    sel_str += "id "+str(all_hydrogens.ids[i])+" or "
sel_str = sel_str[:-4]+")"
sub_hydrogens = u.select_atoms(sel_str)

# logfile per rank
fileout = sys.argv[4]+"_"+str(rank)+".log"

# record current date and time
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

with open(fileout, 'w') as f:
    f.write("Start date and time: "+dt_string+"\n")
    f.write("MPI Size: "+str(size)+"\n")
    f.write("Rank of this process: "+str(rank)+"\n")
    f.write("Atoms in entire universe: "+str(len(u.atoms))+"\n")
    f.write("Residues (molecules) entire in universe: "+str(len(u.residues))+"\n")
    f.write("Number of time steps: "+str(n_ts)+"\n")
    f.write("Atoms per process: "+str(atoms_per_process)+"\n")
    if rank == 0:
        f.write("Sel str: "+sel_str+"\n")
    f.write("Hydrogens to run on this process: "+str(sub_hydrogens.n_atoms)+"\n")
    f.write("Residues to run on this process: "+str(sub_hydrogens.resids)+"\n")
    f.write("Atoms to run on this process: "+str(sub_hydrogens.ids)+"\n")
    f.write("Start and stop range on this process: "+str(startIndex)+" to "+str(endIndex)+"\n")

# run intermolecular autocorrelation calculation
with open(fileout, 'a') as f:
    f.write("Rank "+str(rank)+" START.\n")
inter_acf_local = calc_intra(u, all_hydrogens, sub_hydrogens, n_ts, lag_fraction, fileout) # restore this part for other numpy save file arguments , sys.argv[6]+"_"+str(rank)+".npy", sys.argv[7]+"_"+str(rank)+".npy")

now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
with open(fileout, 'a') as f:
    f.write("Rank "+str(rank)+" DONE.\n")
    f.write("Local ACF: "+str(inter_acf_local)+"\n")
    f.write("Complete date and time: "+dt_string+"\n")
# acf_store_MPI = Comm.gather(inter_acf_local,root=0)

acf_final = Comm.reduce(inter_acf_local, op=MPI.SUM, root=0)/(len(all_hydrogens)/mol_fraction)


if Comm.Get_rank() == 0:
    print("intra_acf:"+str(acf_final))
    with open(fileout, 'a') as f:
        f.write("Final ACF: "+str(acf_final)+"\n")
    np.save(sys.argv[3], acf_final)
    with open(sys.argv[5], 'a') as f:
        f.write('script called: MPI_intra_test.py\n')
        f.write('script run time = '+time.strftime('%H:%M:%S', time.gmtime(time.time()-start))+'\n')
        f.write('acf calculation parameters:\n')
        f.write('fraction of max lag time to total simulation length = '+str(lag_fraction)+'\n')
        f.write('fraction of hydrogens used = 1\n\n')

else:
    pass
