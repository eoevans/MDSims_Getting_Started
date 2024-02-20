import sys
import numpy as np
# import numpy.fft as fft
import MDAnalysis as mda
from datetime import datetime
from mpi4py import MPI
import time
import configparser
from params_validate import *
from os import remove

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
    fvi = np.fft.rfft(vals,n=2*N,axis=0)
    single_acf = np.real(np.fft.irfft(fvi*np.conjugate(fvi),n=2*N,axis=0) [:N,:])
    d = N - np.arange(N)
    d = d[:,np.newaxis]
    d = np.tile(d,(1,single_acf.shape[1]))
    single_acf = single_acf / d 
    single_acf = single_acf[0:int(lagf*len(single_acf))+1,:]
    single_acf = np.sum(single_acf,axis=1)
    return single_acf

def calc_inter(u, all_hydrogens, sub_hydrogen_ids, n_ts, logfile, max_pairs_len):
    # start the looped analysis
    # initialize empty acf_hold for distributed calculation
    acf_local = 0
    # choose a spin
    for i in sub_hydrogen_ids:
        dt_string = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        with open('MDSim/logs/inter_escaped.log', 'a') as f:
            f.write("Hydrogen number "+i+" START at "+dt_string+'\n')
        selH = u.select_atoms('id '+i)
        selH = selH[0]
        # subset other H atoms
        pairs_positions = all_hydrogens.select_atoms("not id "+str(selH.id)+" and not resid "+str(selH.resid)) #intermolecular
        # split up the pairs into smaller sections. This is done to limit memory usage with large simulation cells
        num_pairs_sections = (len(pairs_positions)//max_pairs_len)+1
        indices = np.linspace(0,len(pairs_positions),num_pairs_sections+1).astype(int)
        for i in range(num_pairs_sections):
            section_positions = pairs_positions[indices[i]:indices[i+1]]
            pairs_vals = np.empty([n_ts,len(section_positions)])
            # find 3theta term for each time step
            for ts in u.trajectory:
                pairs_vals[ts.frame,:] = three_theta(selH.position, section_positions.positions)
            acf_local += calc_single_acf(pairs_vals,params['lag_fraction'])
    return acf_local

def write_by_rank(rank, acf):
    np.save(files['acf_log']+'inter_acf'+str(rank)+'.npy', acf)

def combine_delete(size):
    combined_acf = 0
    for i in range(size):
        acf_file = files['acf_log']+'inter_acf'+str(i)+'.npy'
        try:
            combined_acf += np.load(acf_file)
            print('rank '+str(i)+'added to total')
        except:
            pass
    return combined_acf

def get_escaped():
    escaped = []
    with open(files['acf_log']+'inter_MPI_log_0.log', 'r') as f:
        lines = f.readlines()
    for line in lines:
        if line.split(':')[0] == 'MPI Size':
            old_mpi_size = int(line.split(':')[1][1:])
    for i in range(old_mpi_size):
        acf_file = files['acf_log']+'inter_acf'+str(i)+'.npy'
        if os.path.exists(acf_file) == False:
            print(i)
            escaped = escaped + read_escaped(files['acf_log']+'inter_MPI_log_'+str(i)+'.log')
    hydrogens_to_run = ((len(all_hydrogens)*params['mol_fraction'])//old_mpi_size)*old_mpi_size
    return escaped, hydrogens_to_run

def read_escaped(file):
    with open(file, 'r') as f:
        text = f.read()
    return text.split(']')[0].split('[')[1].split()

if __name__ == '__main__':
    with open('MDSim/logs/inter_escaped.log', 'w') as f:
        f.write('step1\n')


    start = time.time()
    files, params = get_files_params(sys.argv[1], 'inter')

    #check_error(files['dcd'])

    np.set_printoptions(threshold=3000)

    alpha = 1 
    ts_fraction = 1 #1/fraction of total time vector to sample
    surface_dims = get_surface_dims(files['surface_details'])
    box_dims = np.array([surface_dims[0]+params['x_len'], surface_dims[1]+params['y_len'], surface_dims[2]+params['z_len']])
    # I ran out of memory when doing the IFFT(FFT). My solution is to chop up the array that is [time_steps, pairs]
    # into smaller arrays, such that not all pairs are calculated at once. Here is the max length of aforementioned array's pairs dimension
    # this is really only a problem with vortex c18a nodes. The current value is for use with 5 processes per c18a node
    max_pairs_len = 3000

    Comm = MPI.COMM_WORLD
    rank = Comm.Get_rank()
    size = Comm.Get_size()

    #load PDB [1] and DCD [2] file for reading in
    u = mda.Universe(files['pdb_sim'],files['dcd'])
    n_ts = len(u.trajectory)//ts_fraction

    #select H atoms on all molecules
    all_hydrogens = u.select_atoms("type H")
    inter_acf_local = np.empty(int(params['lag_fraction']*n_ts)+1)
    acf_final = np.empty(int(params['lag_fraction']*n_ts)+1)

    # Get the total number of hydrogens and divide it among the processors
    # Note: less hydrogens are run than is specified by the user in mol_fraction,
    # as the number of hydrogens to run should be a multiple of the MPI size
    escaped, hydrogens_to_run = get_escaped()
    #print(escaped)
    min_atoms_per_process = len(escaped) // size
    remainder_atoms = len(escaped) % size
    start_index = rank * min_atoms_per_process
    end_index = (rank + 1) * min_atoms_per_process
    if rank >= size - remainder_atoms:
        end_index += rank - (size - remainder_atoms) + 1
        start_index += rank - (size - remainder_atoms)

    # define the atoms that each process will run
    sub_hydrogen_ids = escaped[start_index:end_index]

    # logfile per rank
    fileout = files['acf_log']+'inter_MPI_log_'+str(rank)+'.log'

    # record current date and time
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

    with open('MDSim/logs/inter_escaped.log', 'a') as f:
        f.write('start rank '+str(rank)+'\n')
    inter_acf_local = calc_inter(u, all_hydrogens, sub_hydrogen_ids, n_ts, fileout, max_pairs_len) 
    with open('MDSim/logs/inter_escaped.log', 'a') as f:
        f.write('finished rank '+str(rank)+'\n')

    #np.save(files['acf_log']+'inter_acf'+str(rank)+'.npy', inter_acf_local)

    #Comm.Barrier()
    # for a reason outside my knowledge, this line of code causes the program to hang with large boxes.
    # I might ask the HPC people about this and implement a fix that does not compromise my morals. 
    # But for now, let's have each rank of the process write out a .npy and have rank 0 go back to combine and delete them all.
    acf_final = Comm.reduce(inter_acf_local, op=MPI.SUM, root=0)


    if Comm.Get_rank() == 0:
        np.save(files['acf_log']+'inter_escaped.npy', acf_final)
        print(acf_final[0])
        acf_final = (combine_delete(size)+acf_final)/hydrogens_to_run
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
