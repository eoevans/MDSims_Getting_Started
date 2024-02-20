import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import math
import configparser
from params_validate import *

# Write the taus, dws, and T12s to the results file and graph the acfs. 
# for use with a single trajectory. Updated to handle m=0,1,2 acfs

def graph_acfs(acf_R, acf_T, spacing, file):#    Acf_intra = Acf_intra[:len(t_intra)]
    #Acf_inter = Acf_inter[:len(t_inter)]
    t_axis = np.arange(np.shape(acf_R)[1]) * spacing
    figure, axis = plt.subplots(1,2, figsize=(15,8))
    axis[0].plot(t_axis,acf_R[0,:], t_axis,acf_R[1,:], t_axis,acf_R[2,:])
    axis[0].set_yscale('log')
    axis[0].set_xlabel('t')
    axis[0].set_ylabel('intra_ACF')
    axis[0].legend(['m=0', 'm=1', 'm=2'])
    axis[1].plot(t_axis,acf_T[0,:], t_axis,acf_T[1,:], t_axis,acf_T[2,:])
    axis[1].set_yscale('log')
    axis[1].set_xlabel('t')
    axis[1].set_ylabel('inter_ACF')
    axis[1].legend(['m=0', 'm=1', 'm=2'])
    #axis[0].set_xlim([0,10])
    #axis[0].set_ylim([0.01,1])ax.legend(['First line', 'Second line'])

    plt.savefig(get_sim_dir(file)+'results/acf_graph')
    
def calc_time_scale(spacing, tau):
    # calculate the correct time scale for graphing from t=0 to t=10tau
    return np.arange(0,10*tau,spacing)/tau

def even(arr):
    # make the ACF an even function. This is useful for calculating the spectral density, J, by fourier transform. 
    arr = np.concatenate((np.flip(arr)[0:-1], arr[0:-1]))
    return arr

def calc_J(acf):
    # calculate the spectral density, J, by fourier transform
    return np.fft.rfft(even(acf)).real
    
def Riemann_sum(acf, spacing):
    # calculate the correlation time, tau, by Riemann sum. This is done to elimate confusion regarding the frame spacing
    # when calculating tau with fourier transform
    rsum = np.zeros(3)
    for i in range(np.shape(acf)[1]-1):
        rsum += ((acf[:,i] + acf[:,i+1])/2)*spacing
    return rsum

def carry(hr, min, sec):
    if sec >= 60:
        min += sec // 60
        sec = sec % 60
    if min >= 60:
        hr += min // 60
        min = min % 60
    hr = str(hr)
    min = str(min)
    sec = str(sec)
    if len(hr) == 1:
        hr = '0' + hr
    if len(min) == 1:
        min = '0' + min
    if len(sec) == 1:
        sec = '0' + sec
    return [hr, min, sec]

def values_of_interest(file, numpy_arr=False, return_acf=True):
    if numpy_arr:
        acf = file
    else:
        acf = np.load(file)
    CONSTANT = 1.06805*(10**-13)
    CONSTANTS = np.array([CONSTANT, CONSTANT/(5/(16*np.pi))*(15/(8*np.pi)), CONSTANT/(5/(16*np.pi))*(15/(32*np.pi))])
    framespacing = 0.1
    for m in range(3):
        acf[m,:] = acf[m,:]*CONSTANTS[m]
    tau = Riemann_sum(acf, framespacing)/acf[:,0]
    dw = (((3*acf[:,0])**0.5)*10**9)/(2*np.pi)
    #t12 = ((10*tau*acf[0])**-1)/(10**12)
    T1 = ((2*acf[1,0]*tau[1]+8*acf[2,0]*tau[2])**-1)*(10**-12)
    T2 = ((3*acf[0,0]*tau[0]+5*acf[1,0]*tau[1]+2*acf[2,0]*tau[2])**-1)*(10**-12)
    #graph_acfs(acf, framespacing, file)
    if return_acf:
        return [tau, dw, T1, T2], acf
    else:
        return [tau, dw, T1, T2]

def write_voi(file, acf, numpy_arr=False):
    voi, acf = values_of_interest(acf, numpy_arr=numpy_arr)
    with open(file, 'a') as f:
        f.write('Tau (ps) = ' + str(voi[0]) + '\n')
        f.write('dw (kHz) = '+str(voi[1])+'\n')
        f.write('T1 = ' + str(voi[2]) + '\n')
        f.write('T2 = ' + str(voi[3]) + '\n')

if __name__ == '__main__':
    start = time.time()
    files = read_config(sys.argv[1], 'files')
    # voi is array [[m0tau, m1tau, m2tau], [m0dw, m1dw, m2dw], T1, T2]
    voi_R, acf_R = values_of_interest(files['intra_acf'])
    voi_T, acf_T = values_of_interest(files['inter_acf'])
    framespacing = 0.1
    graph_acfs(acf_R, acf_T, framespacing, files['intra_acf'])

    #u = mda.Universe(,'\\files.campus.wm.edu\acstore-groups\ISC1026\Data\MDSim\data\simulation_results\C5H12\testrun\PNT_box_1.dcd')
    #nframes = len(u.trajectory)
    #framespacing = round(u.trajectory.time, 5)

    total_hr=0
    total_min=0
    total_sec=0
    with open(files['results'], 'r') as f:
        lines = f.readlines()
        for i in lines:
            if i[:15] == 'script run time':
                total_hr += int(i[18:20])
                total_min += int(i[21:23])
                total_sec += int(i[24:26])

    run_time = carry(total_hr, total_min, total_sec)
        
    with open(files['results'], 'a') as f:
        f.write('script called: CalcT12.py\n')
        f.write('script run time = '+time.strftime('%H:%M:%S', time.gmtime(time.time()-start))+'\n')
        f.write('TauR (ps) = ' + str(voi_R[0]) + '\n')
        f.write('TauT (ps) = ' + str(voi_T[0]) + '\n')
        f.write('dwR (kHz) = '+str(voi_R[1])+'\n')
        f.write('dwT (kHz) = '+str(voi_T[1])+'\n')
        f.write('Intramolecular T1 = ' + str(voi_R[2]) + '\n')
        f.write('Intermolecular T1 = ' + str(voi_T[2]) + '\n')
        f.write('Intramolecular T2 = ' + str(voi_R[3]) + '\n')
        f.write('Intermolecular T2 = ' + str(voi_T[3]) + '\n')
        f.write('Total T1 = ' + str(((voi_R[2]**-1)+(voi_T[2])**-1)**-1) + '\n')
        f.write('Total T2 = ' + str(((voi_R[3]**-1)+(voi_T[3])**-1)**-1) + '\n')
        #f.write('Successful graph? ')
        #if graph:
        #    f.write('yes\n\n')
        #else:
        #    f.write('no\n\n')
        f.write('total run time = '+run_time[0]+':'+run_time[1]+':'+run_time[2])