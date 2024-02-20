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
from CalcT12_isotropic import values_of_interest_iso
from CalcT12_anisotropic import values_of_interest
from MDAnalysis import transformations
import matplotlib.pyplot as plt

if __name__ == '__main__':
    files, params = get_files_params(sys.argv[1], 'inter')
    surface_dims = get_surface_dims(files['surface_details'])
    u = mda.Universe(files['pdb_sim'],files['dcd'])
    hydrogens = u.select_atoms('type H')    
    fig, ax = plt.subplots()
    ax.set_xlabel('z position (A)', fontsize=16)
    ax.set_ylabel('probability', fontsize=16)
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='x', labelsize=12)
    # graphing params
    sample_spacing = 0.1
    min = surface_dims[2]-1
    max = params['z_len']+surface_dims[2]+1
    counts, bins = np.histogram(hydrogens.positions[:,2], bins=np.arange(min,max,sample_spacing))
    #bins = bins - bins[0]
    counts = counts/len(hydrogens)
    for ts in u.trajectory[1:]:
        tmp_counts, tmp_bins = np.histogram(hydrogens.positions[:,2], bins=bins)
        counts += tmp_counts/len(hydrogens)
    ax.stairs(counts/len(u.trajectory), bins, fill=True)
    plt.savefig(get_sim_dir(sys.argv[1])+'results/z_pos_prob')
    #counts = counts - counts[len(counts)//2]
    #w = fft.fft(counts)
    #freqs = fft.fftfreq(counts.size)*(sample_spacing**-1)
    #print(freqs[np.argmax(np.abs(w))])
    #fig, ax = plt.subplots()
    #ax.plot(freqs, w)
    #plt.savefig('MDSim/data/numbers_and_graphs/surface_frequency')


