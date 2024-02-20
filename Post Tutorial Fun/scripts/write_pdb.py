import MDAnalysis as mda

dir_sim = 'MDSim/data/simulation_results/C14H30/spread_6-15-23/'
dir_box = 'MDSim/data/box_data/C14H30/'

for i in range(10):
    u = mda.Universe(dir_box+'TDC_box_'+str((i%5)+1)+'.pdb',dir_sim+'trajectory'+str(i+1)+'.dcd')
    u.trajectory[len(u.trajectory)-1]
    u.atoms.write(dir_sim+'run'+str(i+1)+'.inpcrd')