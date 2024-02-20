# Write the results including taus, dws, and T12s, and graph the acfs
# This is currently not updated for use with parameter and file .ini files. 
# Indended for use with a single trajectory. 
# I should be getting rid of this in favor of collect_results.py which puts the reuslts of a single automation run in the same spot.

import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import math
import configparser
from params_validate import *
#import MDAnalysis as mda

start = time.time()

def graph_acfs(Acf_intra, Acf_inter, tau_intra, tau_inter, spacing):#    Acf_intra = Acf_intra[:len(t_intra)]
    #Acf_inter = Acf_inter[:len(t_inter)]
    Acf_intra = Acf_intra / Acf_intra[0]
    Acf_inter = Acf_inter / Acf_inter[0]
    t_intra = np.arange(0,len(Acf_intra))*spacing/tau_intra
    t_inter = np.arange(0,len(Acf_inter))*spacing/tau_inter
    figure, axis = plt.subplots(1,2, figsize=(17,8))
    axis[0].plot(t_intra,Acf_intra)
    axis[0].set_yscale('log')
    axis[0].set_xlabel('t/tr')
    axis[0].set_ylabel('Acf/Acf[0]')
    axis[0].set_title('Intra')
    axis[0].set_xlim([0,10])
    axis[0].set_ylim([0.01,1])
    axis[1].plot(t_inter,Acf_inter)
    axis[1].set_yscale('log')
    axis[1].set_ylabel('Acf/Acf[0]')
    axis[1].set_xlabel('t/tr')
    axis[1].set_title('Inter')
    axis[1].set_xlim([0,10])
    axis[1].set_ylim([0.01,1])
    plt.savefig(files['acf_log'][:-5]+'/results/acf_graph')
    
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
    rsum = 0
    for i in range(len(acf)-1):
        rsum += ((acf[i] + acf[i+1])/2)*spacing
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

def values_of_interest_iso(acf):
    CONSTANT = 1.06805*(10**-13)
    framespacing = 0.1
    print(type(acf))
    acf = acf*CONSTANT
    tau = Riemann_sum(acf, framespacing)/acf[0]
    dw = (((3*acf[0])**0.5)*10**9)/(2*math.pi)
    t12 = ((10*tau*acf[0])**-1)/(10**12)
    return [tau, dw, t12]

if __name__ == '__main__':
    files = read_config(sys.argv[1], 'files')
    Acf_intra = np.load(files['intra_acf'])
    Acf_inter = np.load(files['inter_acf'])

    print(Acf_intra[0])
    print(Acf_inter[0])
    # constant in A^6 ps^-2
    CONSTANT = 1.06805*(10**-13) #this is the constant for m=0

    Acf_intra = Acf_intra*CONSTANT
    Acf_inter = Acf_inter*CONSTANT

    #u = mda.Universe(,'\\files.campus.wm.edu\acstore-groups\ISC1026\Data\MDSim\data\simulation_results\C5H12\testrun\PNT_box_1.dcd')
    #nframes = len(u.trajectory)
    #framespacing = round(u.trajectory.time, 5)
    
    # frame spacing in ps
    framespacing = 0.1

    #Acf_intra = Acf_intra[:1500]
    #Acf_inter = Acf_inter[:1500]
    tau_intra = Riemann_sum(Acf_intra, framespacing)/Acf_intra[0]
    tau_inter = Riemann_sum(Acf_inter, framespacing)/Acf_inter[0]

    print('max lag time intra = '+str(np.size(Acf_intra)*framespacing)+' ps')
    print('size of 10*tau_intra = '+str(tau_intra*10)+' ps')
    print('max lag time inter = '+str(np.size(Acf_inter)*framespacing)+' ps')
    print('size of 10*tau_inter = '+str(tau_inter*10)+' ps')

    dw_R = (((3*Acf_intra[0])**0.5)*10**9)/(2*math.pi)
    dw_T = (((3*Acf_inter[0])**0.5)*10**9)/(2*math.pi)

    T12 = ((10*(tau_intra*Acf_intra[0] + tau_inter*Acf_inter[0]))**-1)/(10**12)
    T12_intra = ((10*tau_intra*Acf_intra[0])**-1)/(10**12)
    T12_inter = ((10*tau_inter*Acf_inter[0])**-1)/(10**12)

    # Check acf decay over entire time
    #np.arange(0,len(Acf_inter),0.1)
    #np.arange(0,len(Acf_intra),0.1)
    #graph = True
    try:
        graph_acfs(Acf_intra,Acf_inter,tau_intra,tau_inter,framespacing)
    except ValueError:
        t_intra = np.arange(0,len(Acf_intra)*framespacing,framespacing)/tau_intra
        t_inter = np.arange(0,len(Acf_inter)*framespacing,framespacing)/tau_inter
        graph_acfs(Acf_intra,Acf_inter,t_intra,t_inter, framespacing)
        print('Error when graphing. Ensure that 10*tau < Max lag time')
        

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
        f.write('Intramolecular Correlation Time (ps) = ' + str(tau_intra) + '\n')
        f.write('Intermolecular Correlation Time (ps) = ' + str(tau_inter) + '\n')
        f.write('dwR (kHz) = '+str(dw_R)+'\n')
        f.write('dwT (kHz) = '+str(dw_T)+'\n')
        f.write('Intramolecular T12 = ' + str(T12_intra) + '\n')
        f.write('Intermolecular T12 = ' + str(T12_inter) + '\n')
        f.write('Total T12 = ' + str(T12) + '\n')
        #f.write('Successful graph? ')
        #if graph:
        #    f.write('yes\n\n')
        #else:
        #    f.write('no\n\n')
        f.write('total run time = '+run_time[0]+':'+run_time[1]+':'+run_time[2])