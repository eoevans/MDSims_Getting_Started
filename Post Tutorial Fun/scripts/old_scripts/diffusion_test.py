# Calculate the diffusion coefficient using equations from Singer. 
import sys
import configparser
from params_validate import validate_section
import MDAnalysis as mda
import numpy as np
#import MDAnalysis.analysis.msd as msd
from MDAnalysis.transformations import nojump
import matplotlib.pyplot as plt
from scipy.stats import linregress
import math

files = configparser.ConfigParser()
files.read(sys.argv[1])
files = files['files']

#u = mda.Universe(files['pdb'],files['dcd'], transformations=nojump.NoJump())
#u = mda.Universe(files['pdb'],files['dcd'])
#atom = u.atoms[100]
#for ts in u.trajectory[::100]:
#    print(atom.position)
#MSD = msd.EinsteinMSD(u, select='all', msd_type='xyz', fft=True)
#MSD.run()
#pls = MSD.results.timeseries
#np.save('MDSim/data/numbers_and_graphs/NVE/msd1.npy', pls)

msd = np.load('MDSim/data/numbers_and_graphs/NVE/msd1.npy')
nframes = len(msd)
timestep = 0.1 # this needs to be the actual time between frames
lagtimes = np.arange(nframes)*timestep # make the lag-time axis
fig = plt.figure()
ax = plt.axes()
# plot the actual MSD
ax.plot(lagtimes, msd)
ax.set_xlabel('lag time (ps)')
ax.set_ylabel('MSD (A^2)')
plt.savefig('MDSim/diffusion')

start_index = 3000
end_index = 18000
linear_model = linregress(lagtimes[start_index:end_index],msd[start_index:end_index])
slope = linear_model.slope
error = linear_model.stderr
# dim_fac is 3 as we computed a 3D msd with 'xyz'
Dsim = slope/6
# Dsim is currently in A^2/ps. Need to convert to m^2/s
Dsim = Dsim*((10**-10)**2)*(10**12)
print(Dsim)
# boltzmann constant in J/K
kb = 1.380649*(10**-23)
# temperature in K
T = 293
# constant (arises also in Madelung constant)
squiggle = 2.837297
# L in meters
L = 40*(10**-10)
# shear viscosity from singer table 2 in cP (centipoise), converted to Pa*s
shear_viscosity = 2.34/1000
DT = Dsim + (kb*T*squiggle/(6*math.pi*shear_viscosity*L))
print(DT)