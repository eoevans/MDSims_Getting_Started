# Collect and tabulate the results of a single automation run. 
# Serves to replace CalcT12.py.
# Now tabulates T12s and diffusion constants

import numpy as np
from tabulate import tabulate
import matplotlib.pyplot as plt
import math
import sys
from MDSim_automate import *
import configparser
from statistics import stdev, mean
from params_validate import *
import traceback
import csv
from datetime import datetime
from CalcT12_anisotropic import values_of_interest

def get_row_of_arr(arr, indeces):
    to_return = []
    for i in range(len(arr)):
        if len(indeces) == 1:
            to_return.append(arr[i][indeces[0]])
        elif len(indeces) == 2:
            to_return.append(arr[i][indeces[0]][indeces[1]])
    return to_return

def make_single_dir(path):
    if os.path.isdir(path) == False:
        os.mkdir(path)
    return path

def Riemann_sum(acf, spacing):
    rsum = 0
    for i in range(len(acf)-1):
        rsum += ((acf[i] + acf[i+1])/2)*spacing
    return rsum    

def str_to_int(to_convert):
    print(to_convert)
    tmp = int(to_convert)
    count = 1
    while tmp == 0:
        print(to_convert[count:])
        tmp = int(to_convert[count:])
        count += 1
    return tmp

def str_to_float(to_convert):
    tmp = to_convert.split('e')
    print(tmp)
    return float(tmp[0]) * (10**str_to_int(tmp[1]))
    
def graph(x,y1,title):
    # y1 is singer, which is short
    figure, axes = plt.subplots(figsize=(10,10))
    axes.plot(x[:len(y1)],y1)
    axes.set_ylabel('Singer/us')
    axes.set_xlabel('Number of Carbons')
    axes.set_title(title)
    plt.savefig('MDSim/data/simulation_results/2ns_6-9-23_results/'+title)

def check_singer(mols, table_type):
    no_singer = {table_type: mols, 'T1 (s)': np.zeros(len(mols)),'T2 (s)': np.zeros(len(mols)),
                'T1,R (s)': np.zeros(len(mols)), 'T1,T (s)': np.zeros(len(mols)),
                'T2,R (s)': np.zeros(len(mols)), 'T2,T (s)': np.zeros(len(mols)),
                'tauR m=0 (ps)': np.zeros(len(mols)),'tauT m=0 (ps)': np.zeros(len(mols)),
                'tauR m=1 (ps)': np.zeros(len(mols)),'tauT m=1 (ps)': np.zeros(len(mols)),
                'tauR m=2 (ps)': np.zeros(len(mols)),'tauT m=2 (ps)': np.zeros(len(mols)),
                'omegaR m=0 (kHz)': np.zeros(len(mols)),'omegaT m=0 (kHz)': np.zeros(len(mols)),
                'omegaR m=1 (kHz)': np.zeros(len(mols)),'omegaT m=1 (kHz)': np.zeros(len(mols)),
                'omegaR m=2 (kHz)': np.zeros(len(mols)),'omegaT m=2 (kHz)': np.zeros(len(mols))}
    return no_singer

def make_T12_table(folders, rows, table_type, table_folder, table_name):
    if table_type == 'mols':
        table_type = 'molecule'
    elif table_type == 'params':
        table_type = 'parameters'
    table_indiv = {table_type: ['']*(len(rows)*len(folders[0])+len(rows)),'T1 (s)': np.zeros(len(rows)*len(folders[0])+len(rows)),
                   'T2 (s)': np.zeros(len(rows)*len(folders[0])+len(rows)),
                   'T1,R (s)': np.zeros(len(rows)*len(folders[0])+len(rows)), 'T1,T (s)': np.zeros(len(rows)*len(folders[0])+len(rows)),
                   'T2,R (s)': np.zeros(len(rows)*len(folders[0])+len(rows)), 'T2,T (s)': np.zeros(len(rows)*len(folders[0])+len(rows)),
                   'tauR m=0 (ps)': np.zeros(len(rows)*len(folders[0])+len(rows)),'tauT m=0 (ps)': np.zeros(len(rows)*len(folders[0])+len(rows)),
                   'tauR m=1 (ps)': np.zeros(len(rows)*len(folders[0])+len(rows)),'tauT m=1 (ps)': np.zeros(len(rows)*len(folders[0])+len(rows)),
                   'tauR m=2 (ps)': np.zeros(len(rows)*len(folders[0])+len(rows)),'tauT m=2 (ps)': np.zeros(len(rows)*len(folders[0])+len(rows)),
                   'omegaR m=0 (kHz)': np.zeros(len(rows)*len(folders[0])+len(rows)),'omegaT m=0 (kHz)': np.zeros(len(rows)*len(folders[0])+len(rows)),
                   'omegaR m=1 (kHz)': np.zeros(len(rows)*len(folders[0])+len(rows)),'omegaT m=1 (kHz)': np.zeros(len(rows)*len(folders[0])+len(rows)),
                   'omegaR m=2 (kHz)': np.zeros(len(rows)*len(folders[0])+len(rows)),'omegaT m=2 (kHz)': np.zeros(len(rows)*len(folders[0])+len(rows))}
    table_sum = check_singer(rows, table_type)
    #dir_list_outer = [r'data\simulation_results\C20H42']
    count = 0
    for i in range(len(folders)):
        table_indiv[table_type][count] = rows[i]
        voi_R = [0]*len(folders[0])
        voi_T = [0]*len(folders[0])
        T2 = np.zeros(len(folders[0]))
        T1 = np.zeros(len(folders[0]))
        for j in range(len(folders[i])):
            print(folders[i])
            voi_R[j] = values_of_interest(folders[i][j]+'intra_acf_MPI.npy', numpy_arr=False, return_acf=False)
            voi_T[j] = values_of_interest(folders[i][j]+'inter_acf_MPI.npy', numpy_arr=False, return_acf=False)
            T1[j] = (voi_R[j][2]**-1+voi_T[j][2]**-1)**-1
            T2[j] = (voi_R[j][3]**-1+voi_T[j][3]**-1)**-1

            table_indiv['T1 (s)'][count] = T1[j]
            table_indiv['T2 (s)'][count] = T2[j]
            table_indiv['T1,R (s)'][count] = voi_R[j][2]
            table_indiv['T1,T (s)'][count] = voi_T[j][2]
            table_indiv['T2,R (s)'][count] = voi_R[j][3]
            table_indiv['T2,T (s)'][count] = voi_T[j][3]
            table_indiv['tauR m=0 (ps)'][count] = voi_R[j][0][0]
            table_indiv['tauT m=0 (ps)'][count] = voi_T[j][0][0]
            table_indiv['tauR m=1 (ps)'][count] = voi_R[j][0][1]
            table_indiv['tauT m=1 (ps)'][count] = voi_T[j][0][1]
            table_indiv['tauR m=2 (ps)'][count] = voi_R[j][0][2]
            table_indiv['tauT m=2 (ps)'][count] = voi_T[j][0][2]
            table_indiv['omegaR m=0 (kHz)'][count] = voi_R[j][1][0]
            table_indiv['omegaT m=0 (kHz)'][count] = voi_T[j][1][0]
            table_indiv['omegaR m=1 (kHz)'][count] = voi_R[j][1][1]
            table_indiv['omegaT m=1 (kHz)'][count] = voi_T[j][1][1]
            table_indiv['omegaR m=2 (kHz)'][count] = voi_R[j][1][2]
            table_indiv['omegaT m=2 (kHz)'][count] = voi_T[j][1][2]
            count += 1
        
        csv_dir = make_single_dir('MDSim/data/numbers_and_graphs/'+table_folder+'/'+table_name+'_csvs/')
        with open(csv_dir+'T1.csv', 'a') as f:
            writer = csv.writer(f)
            writer.writerow(T1) 
        with open(csv_dir+'T2.csv', 'a') as f:
            writer = csv.writer(f)
            writer.writerow(T2) 
        with open(csv_dir+'T1_R.csv', 'a') as f:
            writer = csv.writer(f)
            writer.writerow(get_row_of_arr(voi_R,[2]))
        with open(csv_dir+'T1_T.csv', 'a') as f:
            writer = csv.writer(f)
            writer.writerow(get_row_of_arr(voi_T,[2]))
        with open(csv_dir+'T2_R.csv', 'a') as f:
            writer = csv.writer(f)
            writer.writerow(get_row_of_arr(voi_R,[3]))
        with open(csv_dir+'T2_T.csv', 'a') as f:
            writer = csv.writer(f)
            writer.writerow(get_row_of_arr(voi_T,[3]))

        table_indiv[table_type][count] = 'stdev'
        if len(folders[0]) > 1:
            table_indiv['T1 (s)'][count] = stdev(T1)
            table_indiv['T2 (s)'][count] = stdev(T2)
            table_indiv['T1,R (s)'][count] = stdev(get_row_of_arr(voi_R,[2]))
            table_indiv['T1,T (s)'][count] = stdev(get_row_of_arr(voi_T,[2]))
            table_indiv['T2,R (s)'][count] = stdev(get_row_of_arr(voi_R,[3]))
            table_indiv['T2,T (s)'][count] = stdev(get_row_of_arr(voi_T,[3]))
            table_indiv['tauR m=0 (ps)'][count] = stdev(get_row_of_arr(voi_R,[0,0]))
            table_indiv['tauT m=0 (ps)'][count] = stdev(get_row_of_arr(voi_T,[0,0]))
            table_indiv['tauR m=1 (ps)'][count] = stdev(get_row_of_arr(voi_R,[0,1]))
            table_indiv['tauT m=1 (ps)'][count] = stdev(get_row_of_arr(voi_T,[0,1]))
            table_indiv['tauR m=2 (ps)'][count] = stdev(get_row_of_arr(voi_R,[0,2]))
            table_indiv['tauT m=2 (ps)'][count] = stdev(get_row_of_arr(voi_T,[0,2]))
            table_indiv['omegaR m=0 (kHz)'][count] = stdev(get_row_of_arr(voi_R,[1,0]))
            table_indiv['omegaT m=0 (kHz)'][count] = stdev(get_row_of_arr(voi_T,[1,0]))
            table_indiv['omegaR m=1 (kHz)'][count] = stdev(get_row_of_arr(voi_R,[1,1]))
            table_indiv['omegaT m=1 (kHz)'][count] = stdev(get_row_of_arr(voi_T,[1,1]))
            table_indiv['omegaR m=2 (kHz)'][count] = stdev(get_row_of_arr(voi_R,[1,2]))
            table_indiv['omegaT m=2 (kHz)'][count] = stdev(get_row_of_arr(voi_T,[1,2]))
            count += 1

        ind = table_sum[table_type].index(rows[i]) 
        table_sum['T1 (s)'][ind] = mean(T1)
        table_sum['T2 (s)'][ind] = mean(T2)
        table_sum['T1,R (s)'][ind] = mean(get_row_of_arr(voi_R,[2]))
        table_sum['T1,T (s)'][ind] = mean(get_row_of_arr(voi_T,[2]))
        table_sum['T2,R (s)'][ind] = mean(get_row_of_arr(voi_R,[3]))
        table_sum['T2,T (s)'][ind] = mean(get_row_of_arr(voi_T,[3]))
        table_sum['tauR m=0 (ps)'][ind] = mean(get_row_of_arr(voi_R,[0,0]))
        table_sum['tauT m=0 (ps)'][ind] = mean(get_row_of_arr(voi_T,[0,0]))
        table_sum['tauR m=1 (ps)'][ind] = mean(get_row_of_arr(voi_R,[0,1]))
        table_sum['tauT m=1 (ps)'][ind] = mean(get_row_of_arr(voi_T,[0,1]))
        table_sum['tauR m=2 (ps)'][ind] = mean(get_row_of_arr(voi_R,[0,2]))
        table_sum['tauT m=2 (ps)'][ind] = mean(get_row_of_arr(voi_T,[0,2]))
        table_sum['omegaR m=0 (kHz)'][ind] = mean(get_row_of_arr(voi_R,[1,0]))
        table_sum['omegaT m=0 (kHz)'][ind] = mean(get_row_of_arr(voi_T,[1,0]))
        table_sum['omegaR m=1 (kHz)'][ind] = mean(get_row_of_arr(voi_R,[1,1]))
        table_sum['omegaT m=1 (kHz)'][ind] = mean(get_row_of_arr(voi_T,[1,1]))
        table_sum['omegaR m=2 (kHz)'][ind] = mean(get_row_of_arr(voi_R,[1,2]))
        table_sum['omegaT m=2 (kHz)'][ind] = mean(get_row_of_arr(voi_T,[1,2]))

    #os.mkdir('MDSim/data/simulation_results/2ns_6-9-23_results/')
    with open('MDSim/data/numbers_and_graphs/'+table_folder+'/'+table_name+'_T12_avg.txt', 'w') as f:
        f.write(tabulate(table_sum, headers = 'keys', numalign='center', stralign='center'))
        dt_string = datetime.now().strftime('%d/%m/%Y %H:%M:%S')
        f.write('\ncreated '+dt_string)
    with open('MDSim/data/numbers_and_graphs/'+table_folder+'/'+table_name+'_T12_indiv.txt', 'w') as f:
        f.write(tabulate(table_indiv, headers = 'keys', numalign='center', stralign='center'))
        f.write('\ncreated '+dt_string)

def make_diffusion_table(folders, rows, table_type, table_folder, table_name):
    if table_type == 'mols':
        table_type = 'molecule'
    elif table_type == 'params':
        table_type = 'parameters'
    table_indiv = {table_type: ['']*(len(rows)*len(folders[0])+len(rows)),'Dsim (m^2 s^-1)': np.zeros(len(rows)*len(folders[0])+len(rows)),
                   'DT (m^2 s^-1)': np.zeros(len(rows)*len(folders[0])+len(rows))}
    table_sum = {table_type: rows, 'Dsim (m^2 s^-1)': np.zeros(len(rows)),  'DT (m^2 s^-1)': np.zeros(len(rows))}
    #dir_list_outer = [r'data\simulation_results\C20H42']
    count = 0
    for i in range(len(folders)):
        table_indiv[table_type][count] = rows[i]
        Dsim = np.zeros(len(folders[0]))
        DT = np.zeros(len(folders[0]))
        for j in range(len(folders[i])):
            files = read_config(folders[i][j]+'files.ini', 'files')
            with open(files['results'], 'r') as f:
                lines = f.readlines()
            for line in lines:
                tmp = line.split(' = ')
                if tmp[0] == 'Dsim':
                    Dsim[j] = str_to_float(tmp[1][:-1])
                if tmp[0] == 'DT':
                    DT[j] = str_to_float(tmp[1][:-1])

            table_indiv['Dsim (m^2 s^-1)'][count] = Dsim[j]
            table_indiv['DT (m^2 s^-1)'][count] = DT[j]
            count += 1

        table_indiv[table_type][count] = 'stdev'
        if len(folders[0]) > 1:
            table_indiv['Dsim (m^2 s^-1)'][count] = stdev(Dsim)
            table_indiv['DT (m^2 s^-1)'][count] = stdev(DT)
            count += 1

        ind = table_sum[table_type].index(rows[i]) 
        table_sum['Dsim (m^2 s^-1)'][ind] = mean(Dsim)
        table_sum['DT (m^2 s^-1)'][ind] = mean(DT)

    #os.mkdir('MDSim/data/simulation_results/2ns_6-9-23_results/')
    with open('MDSim/data/numbers_and_graphs/'+table_folder+'/'+table_name+'_diffusion_avg.txt', 'w') as f:
        f.write(tabulate(table_sum, headers = 'keys', numalign='center', stralign='center'))
    with open('MDSim/data/numbers_and_graphs/'+table_folder+'/'+table_name+'_diffusion_indiv.txt', 'w') as f:
        f.write(tabulate(table_indiv, headers = 'keys', numalign='center', stralign='center'))

def change_dimensions(folders):
    params_dim = [None]*len(folders[0])
    for i in range(len(folders[0])):
        mol_dim = [None]*len(folders)
        for j in range(len(folders)):
            run_dim = [None]*len(folders[0][0])
            for k in range(len(folders[j][i])):
                run_dim[k] = folders[j][i][k]
            mol_dim[j] = run_dim
        params_dim[i] = mol_dim
    return params_dim

if __name__ == '__main__':
    # some constants
    CONSTANT = 1.06805*(10**-13)
    framespacing = 0.1
    # read in the necessary sections of the automation parameters
    automate = read_config(sys.argv[1]+'auto_params.ini', 'automate')
    results = read_config(sys.argv[1]+'auto_params.ini', 'results_sub')
    # read in the molecule array
    mols = get_mols(automate)
    params_folders = get_params_folders(sys.argv[1])
    # get the simulation folders, i.e. where the trajectory and results are stored. 
    sim_folders = get_folders(automate,mols,params_folders[1],'sim_folder')
    if results['table_rows'] == 'mols':
        # if mols is specified as a rows, generate n tables, where n is the number of parameters. 
        sim_folders = change_dimensions(sim_folders)
        for i in range(len(sim_folders)):
            try:
                make_T12_table(sim_folders[i], mols, results['table_rows'], automate['sim_folder'], params_folders[1][i])
            except:
                print('error in T12 table')
                traceback.print_exc()
            try:
                make_diffusion_table(sim_folders[i], mols, results['table_rows'], automate['sim_folder'], params_folders[1][i])
            except:
                print('error in diffusion table')
                traceback.print_exc()
    elif results['table_rows'] == 'params':
        # if params is specified as a rows, generate n tables, where n is the number of molecules.
        for i in range(len(sim_folders)):
            try:
                make_T12_table(sim_folders[i], params_folders[1], results['table_rows'], automate['sim_folder'], mols[i])
            except:
                print('error in T12 table')
                traceback.print_exc()
            try:
                make_diffusion_table(sim_folders[i], params_folders[1], results['table_rows'], automate['sim_folder'], mols[i])
            except:
                print('error in diffusion table')
                traceback.print_exc()
    
    
    

    #Singer_len = len(table['T12 (s) Singer'])
    #x_axis = [5,6,8,10,12,14,16,18,20]
    #graph(x_axis,table['T12 (s) Singer']/table['T12 (s) us'][:Singer_len], 'T12')
    #graph(x_axis,table['T12,R (s) Singer']/table['T12,R (s) us'][:Singer_len], 'T12,R')
    #graph(x_axis,table['T12,T (s) Singer']/table['T12,T (s) us'][:Singer_len], 'T12')
    #graph(x_axis,table['tauR (ps) Singer']/table['tauR (ps) us'][:Singer_len], 'tauR')
    #graph(x_axis,table['tauT (ps) Singer']/table['tauT (ps) us'][:Singer_len], 'tauT')
    #graph(x_axis,table['omegaR (kHz) Singer']/table['omegaR (kHz) us'][:Singer_len], 'omegaR')
    #graph(x_axis,table['omegaT (kHz) Singer']/table['omegaT (kHz) us'][:Singer_len], 'omegaT')