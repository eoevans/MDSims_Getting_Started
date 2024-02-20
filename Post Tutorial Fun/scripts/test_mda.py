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

def get_cutoff_distance(box_dims):
    # args are box dimensions in A. Returns half the smallest dimension of the box or 20 A per Singer
    cutoff_distance = 10
    for i in box_dims:
        if i/2 < cutoff_distance:
            cutoff_distance = i/2
    return cutoff_distance

def three_theta(point, points):
    v = np.absolute(points-point)
    v = v - np.trunc(v/box_dims+0.5)*box_dims
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

if __name__ == '__main__':
    files, params = get_files_params(sys.argv[1], 'inter')
    surface_dims = get_surface_dims(files['surface_details'])
    box_dims = np.array([surface_dims[0]+params['x_len'], surface_dims[1]+params['y_len'], surface_dims[2]+params['z_len']])
    u = mda.Universe(files['pdb_sim'],files['dcd'])
    all_hydrogens = u.select_atoms("name H*")
    print(len(all_hydrogens))
    selH = all_hydrogens.select_atoms('index 2891')
    radius = get_cutoff_distance(box_dims)
    print(len(selH))
    for i in range(len(selH)):
        print(selH[i])
        pairs_positions = get_close_hydrogens(all_hydrogens, selH[i], radius)
        print(len(pairs_positions))