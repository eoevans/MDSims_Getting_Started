import numpy as np
import matplotlib.pyplot as plt
import math

def Riemann_sum(acf, spacing):
    # calculate the correlation time, tau, by Riemann sum. This is done to elimate confusion regarding the frame spacing
    # when calculating tau with fourier transform
    rsum = 0
    for i in range(len(acf)-1):
        rsum += ((acf[i] + acf[i+1])/2)*spacing
    return rsum

def calc_time_scale(spacing, tau):
    # calculate the correct time scale for graphing from t=0 to t=10tau
    return np.arange(0,10*tau,spacing)/tau

if __name__ == '__main__':
    Acf_inter = 0
    Acf = 0
    Acf_short = 0
    Acf_T_short = 0
    CONSTANT = 1.06805*(10**-13)
    framespacing = 0.1
    mol = 'C14H30'

    base_dir = 'MDSim/data/simulation_results/'+mol+'/NVE/run'

    type_acf = ['intra','inter']

    for j in type_acf:
        figure, axis = plt.subplots(1,2,figsize=(15,8))
        #legend = ['run1','run2','run3','run4','run5']

        count = -1
        for i in range(5):
            count += 1
            Acf = np.load(base_dir+str(i+1)+'/'+j+'_acf_MPI.npy')

            Acf = Acf*CONSTANT

    #u = mda.Universe(,'\\files.campus.wm.edu\acstore-groups\ISC1026\Data\MDSim\data\simulation_results\C5H12\testrun\PNT_box_1.dcd')
    #nframes = len(u.trajectory)
    #framespacing = round(u.trajectory.time, 5)

    #Acf = Acf[:1500]
    #Acf_inter = Acf_inter[:1500]
            tau = Riemann_sum(Acf, framespacing)/Acf[0]
            dw_R = (((3*Acf[0])**0.5)*10**9)/(2*math.pi)
            T12_R = ((10*tau*Acf[0])**-1)/(10**12)
            print(tau)
        #print(tau)
        #tau_inter = Riemann_sum(Acf_inter, framespacing)/Acf_inter[0]

        #Acf_inter = Acf_inter[:len(t_inter)]
            Acf = Acf / Acf[0]
        #Acf_inter = Acf_inter / Acf_inter[0]
            t = calc_time_scale(framespacing,tau)
        #t_inter = calc_time_scale(framespacing,tau_inter)
            Acf_short = Acf
        #Acf_T_short = Acf_inter
            if 10*tau <= len(Acf)*framespacing:
                Acf_short = Acf[:len(t)]
            else:
                t = t[0:len(Acf)]
        #if 10*tau_inter <= len(Acf_inter)*framespacing:
        #    Acf_T_short = Acf_inter[:len(t_inter)]
        #else:
        #    t_inter = t_inter[0:len(Acf_inter)]

            axis[0].plot(t,Acf_short)
        #axis[1][0].plot(t_inter,Acf_T_short,label=legend[count])
            lagt = np.arange(0,len(Acf)*framespacing,framespacing)
        #lagt_T = np.arange(0,len(Acf_inter)*framespacing,framespacing)
            axis[1].plot(lagt,Acf)
        #axis[1][1].plot(lagt_T,Acf_inter,label=legend[count])
        axis[0].set_yscale('log')
        axis[0].set_xlabel('t/tau_'+j)
        axis[0].set_ylabel('Acf/Acf[0]')
        axis[0].set_title(j)
        axis[0].set_xlim([0,10])
        axis[0].set_ylim([0.01,1])
        axis[1].set_yscale('log')
        axis[1].set_xlabel('lag time (ps)')
        axis[1].set_ylabel('Acf/Acf[0]')
        axis[1].set_title(j)
        axis[1].set_ylim([0.01,1])
        plt.savefig('MDSim/data/numbers_and_graphs/NVE/'+j)
   
    