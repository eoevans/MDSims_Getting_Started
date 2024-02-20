import sys
import numpy as np
# import numpy.fft as fft
import pyfftw.interfaces.numpy_fft as fft
import MDAnalysis as mda
from datetime import datetime
from mpi4py import MPI
import time
import configparser
from params_validate import get_files_params

# make function to find the 3theta term
def three_theta(point, points):
    v = np.absolute(points-point)
    v = v - np.trunc(v/params['side_len']+0.5)*params['side_len']
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

def calc_inter(u, all_hydrogens, sub_hydrogens, n_ts, logfile):
    # how many molecular pairs will each H atom have
    n_pairs = all_hydrogens.n_atoms-len(all_hydrogens.select_atoms("resid 1")) #for intermolecular pairs
    
    # initialize empty array that is (time steps, connections)
    pairs_vals = np.empty([n_ts,n_pairs])

    # start the looped analysis
    # initialize empty acf_hold for distributed calculation
    acf_local = 0
    acf_hold = 0
    # choose a spin
    for i in np.arange(len(sub_hydrogens)):
        selH = sub_hydrogens[i]
        dt_string = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        with open(logfile, 'a') as f:
            f.write("Hydrogen number "+str(sub_hydrogens[i].id)+" START at "+dt_string+".\n")

        # subset other H atoms
        pairs_positions = all_hydrogens.select_atoms("not index "+str(selH.id)+" and not resid "+str(selH.resid)) #intermolecular
        
        # find 3theta term for each time step
        for ts in u.trajectory[0:n_ts-1]:
            pairs_vals[ts.frame,:] = three_theta(selH.position, pairs_positions.positions)
        acf_local += calc_single_acf(pairs_vals,params['lag_fraction'])
    
    #acf_hold = acf_local/n_pairs
    return acf_local

print('starting\n')
start = time.time()
files, params = get_files_params(sys.argv[1], 'inter')

#check_error(files['dcd'])

np.set_printoptions(threshold=3000)

alpha = 1 #2.025568e-15 #m^6 s^-2, this is the constant in front of the whole expression for protons
# alternatively, 2.025568e-15 #A^6 ps^-2
ts_fraction = 1 #1/fraction of total time vector to sample

Comm = MPI.COMM_WORLD
rank = Comm.Get_rank()
size = Comm.Get_size()

#load PDB [1] and DCD [2] file for reading in
u = mda.Universe(files['pdb_box'],files['dcd'])
n_ts = len(u.trajectory)//ts_fraction

#select H atoms on all molecules
all_hydrogens = u.select_atoms("name H*")
inter_acf_local = np.empty(int(params['lag_fraction']*n_ts)+1)
acf_final = np.empty(int(params['lag_fraction']*n_ts)+1)


if rank == 0:
    #select H atoms on a subset the molecules
    #randomly choose the subset of the molecules
    n_choose = len(u.residues)
    n_H = len(u.select_atoms("resid 1").ids)
    while (n_choose*n_H) % size != 0:
        print(n_choose)
        n_choose += 1

    # pick_atom = 2 #use this to pick a spefici atom
    # sub_hydrogens = u.select_atoms("id "+str(pick_atom))
    rand_sel = np.random.permutation(n_choose)
    rand_sel = np.sort(rand_sel[0:int(n_choose*params['mol_fraction'])])
    print("rand_sel:"+str(rand_sel))

    # make string for AtomGroup selection of the subset of the molecules
    sel_str = '('
    for i in np.arange(len(rand_sel)):
        sel_str += "resid "+str(rand_sel[i])+" or "
    sel_str = sel_str[:-4]+") and name H*"
    sub_hydrogen_ids = u.select_atoms(sel_str).ids
    print(f"sub_H ids: {sub_hydrogen_ids}")
else:
    sub_hydrogen_ids = None

# broadcast list of sub-sampled H atoms to whole group
sub_hydrogen_ids = Comm.bcast(sub_hydrogen_ids, root=0)

# Get the total number of hydrogens and divide it among the processors
atoms_per_process = len(sub_hydrogen_ids) // size
rank = Comm.Get_rank()
startIndex = rank * atoms_per_process
endIndex = (rank + 1) * atoms_per_process

# make a selection string for each process to find its own H atoms
sel_str = '('
for i in range(startIndex,endIndex):
    sel_str += "id "+str(sub_hydrogen_ids[i])+" or "
sel_str = sel_str[:-4]+")"

sub_hydrogens = u.select_atoms(sel_str)

# logfile per rank
fileout = files['acf_log']+'inter_MPI_log_'+str(rank)+'.log'

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
inter_acf_local = calc_inter(u, all_hydrogens, sub_hydrogens, n_ts, fileout) # restore this part for other numpy save file arguments , sys.argv[6]+"_"+str(rank)+".npy", sys.argv[7]+"_"+str(rank)+".npy")

now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
with open(fileout, 'a') as f:
    f.write("Rank "+str(rank)+" DONE.\n")
    f.write("Local ACF: "+str(inter_acf_local)+"\n")
    f.write("Complete date and time: "+dt_string+"\n")
# acf_store_MPI = Comm.gather(inter_acf_local,root=0)

acf_final = Comm.reduce(inter_acf_local, op=MPI.SUM, root=0)/(len(all_hydrogens)*params['mol_fraction'])


if Comm.Get_rank() == 0:
    print("inter_acf:"+str(acf_final))
    with open(fileout, 'a') as f:
        # f.write("Gathered ACF: "+str(acf_store_MPI)+"\n")
        f.write("Final ACF: "+str(acf_final)+"\n")
    np.save(files['inter_acf'], acf_final)

    with open(files['results'], 'a') as f:
        f.write('inter\n')
        f.write('script run time = '+time.strftime('%H:%M:%S', time.gmtime(time.time()-start))+'\n')
        for key in params:
            f.write(key+' = '+str(params[key])+'\n')
        f.write('\n')
        
else:
    pass
