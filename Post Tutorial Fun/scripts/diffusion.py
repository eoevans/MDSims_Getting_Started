# Calculate the diffusion coefficient using equations from Singer. 
import sys
import configparser
from params_validate import get_files_params
import MDAnalysis as mda
import numpy as np
import MDAnalysis.analysis.msd as msd
from MDAnalysis.transformations import nojump
import matplotlib.pyplot as plt
from scipy.stats import linregress
import constants
import time

start = time.time()
files, params = get_files_params(sys.argv[1], 'diffusion')

u = mda.Universe(files['pdb_sim'],files['dcd'], transformations=nojump.NoJump())
MSD = msd.EinsteinMSD(u, select='all', msd_type='xyz', fft=True)
MSD.run()
msd_calced = MSD.results.timeseries
#np.save('MDSim/data/numbers_and_graphs/NVE/msd1.npy', pls)

#msd = np.load('MDSim/data/numbers_and_graphs/NVE/msd1.npy')
nframes = len(msd_calced)
timestep = 0.1 # this needs to be the actual time between frames
lagtimes = np.arange(nframes)*timestep # make the lag-time axis
fig = plt.figure()
ax = plt.axes()
# plot the actual MSD
ax.plot(lagtimes, msd_calced)
ax.set_xlim([lagtimes[int(params['drop_start']*nframes)], lagtimes[int(params['drop_end']*nframes)*-1]])
ax.set_xlabel('lag time (ps)')
ax.set_ylabel('MSD (A^2)')
plt.savefig(files['diffusion_graph'])

linear_model = linregress(lagtimes[int(params['drop_start']*nframes):int(params['drop_end']*nframes)*-1],msd_calced[int(params['drop_start']*nframes):int(params['drop_end']*nframes*-1)])
slope = linear_model.slope
error = linear_model.stderr

# Calculate Dsim convert it from A^2/ps to m^2/s
Dsim = (slope/6)*((10**-10)**2)*(10**12)
# assuming simulation cell is a cube
DT = Dsim + (constants.kb*params['temp']*constants.squiggle/(6*constants.pi*constants.shear_vis(files['mol'])*params['x_len']*10**-10))
print(files['results'])
with open(files['results'],'a') as f:
    f.write('diffusion:\n')
    f.write('script run time = '+time.strftime('%H:%M:%S', time.gmtime(time.time()-start))+'\n')
    f.write('Dsim = '+str(Dsim)+'\n')
    f.write('DT = '+str(DT)+'\n\n')