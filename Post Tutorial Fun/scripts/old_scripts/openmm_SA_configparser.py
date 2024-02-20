# This script was generated by OpenMM-Setup on 2023-02-28.

import sys
from openmm import *
from openmm.app import *
from simtk.unit import *
import time
import configparser
from params_validate import validate_section

# start the clock
start = time.time()

# read in parameters
params = validate_section(sys.argv[1], 'openmm')

files = configparser.ConfigParser()
files.read(sys.argv[2])
files = files['files']

# Input Files

prmtop = AmberPrmtopFile(files['prmtop'])
inpcrd = AmberInpcrdFile(files['inpcrd'])

# System Configuration

nonbondedMethod = PME
nonbondedCutoff = 1.4*nanometers
#nonbondedCutoff = 0.7*nanometers
ewaldErrorTolerance = 0.0005
constraints = None
rigidWater = True

# Integration Options

dt = 0.001*picoseconds
friction = 1.0/picosecond
pressure = 1.0*atmospheres
barostatInterval = 25

# Simulation Options
params['temp_step_size'] = (params['initial_temp'] - params['final_temp'])/params['temp_steps']
params['steps_per_temp'] = int(params['equil_steps']/params['temp_steps'])

# intdroduce boolean to turn GPU on or off
if params['gpu']:
    platform = Platform.getPlatformByName('CUDA') #For CUDA, use 'CUDA'
    platformProperties = {'Precision': 'double'} #For CUDA, use 'Precision': 'double' or 'mixed' (for accuracy and speed, respectively.
else:
    platform = Platform.getPlatformByName('CPU')
    platformProperties = {}
# pdbReporter = PDBReporter(sys.argv[5], steps)
dcd_interval = 100
dcdReporter = DCDReporter(files['dcd'], dcd_interval)
dataReporter = StateDataReporter(files['openmm_log'], 1000, totalSteps=params['sim_steps'],
    step=True, time=True, speed=True, progress=True, elapsedTime=True, remainingTime=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, separator='\t')

# Prepare the Simulation

print('Building system...')
topology = prmtop.topology
positions = inpcrd.positions
system = prmtop.createSystem(nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff,
    constraints=constraints, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance)
# system.addForce(MonteCarloBarostat(pressure, temperature, barostatInterval))
integrator = LangevinMiddleIntegrator(params['initial_temp']*kelvin, friction, dt)
# integrator = VerletIntegrator(dt)
simulation = Simulation(topology, system, integrator)
simulation.context.setPositions(positions)
# print(f"inpcrd.boxVectors is {inpcrd.boxVectors}\n")
if inpcrd.boxVectors is not None:
    print(f"inpcrd.boxVectors is {inpcrd.boxVectors}\n")
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

# Minimize and Equilibrate

print('Performing energy minimization...')
simulation.minimizeEnergy()
# print('Equilibrating...')
for i in range(params['temp_steps']):
    integrator.setTemperature((params['initial_temp']-params['temp_step_size']*i)*kelvin)
    simulation.step(params['steps_per_temp'])
# simulation.context.setVelocitiesToTemperature(temperature)
# simulation.step(equilibrationSteps)

# Simulate
print('Simulating...')
integrator.setTemperature(params['final_temp']*kelvin)
simulation.reporters.append(dcdReporter)
simulation.reporters.append(dataReporter)
# simulation.reporters.append(pdbReporter)
simulation.currentStep = 0
simulation.step(params['sim_steps'])
# for i in range(200):
#     integrator.setTemperature(3*(300-i)*kelvin)
#     simulation.step(1000)

# end simulation section -- now to write the simulation parameters to the results file
params['num_trajectory_frames'] = params['sim_steps']/dcd_interval
params['run_time'] = params['sim_steps']*dt
params['frame_spacing'] = params['run_time']/params['num_trajectory_frames']
params['equil_time'] = params['equil_steps']*dt

with open(files['results'], 'w') as f:
    f.write('Openmm\n')
    f.write('All temperatures are in Kelvin.\n')
    #f.write('script called: openmm_sim_anneal.py\n')
    #f.write('script run time = '+time.strftime('%H:%M:%S', time.gmtime(time.time()-start))+'\n')
    #f.write('GPU used? ')
    #if params['GPU']:
    #    f.write('yes\n')
    #else:
    #    f.write('no\n')
    #f.write('simmulated annealing parameters: \n')
    #f.write('start temp: '+str(initial_temperature)+'\n')
    #f.write('end temp: '+str(final_temp)+'\n')
    #f.write('equilibration time: '+str(equil_time)+'\n')
    #f.write('temperature step size: '+str(temp_step_size)+'\n')
    #f.write('time at each temperature: '+str(num_steps_per_temp*dt)+'\n')
    #f.write('simulation parameters: \n')
    #f.write('time step = '+str(dt)+'\n')
    #f.write('simulation length = '+str(run_time)+'\n')
    #f.write('number of frames = '+str(int(num_frames))+'\n')
    #f.write('frame spacing = '+str(frame_spacing)+'\n\n')
    for key in params:
        f.write(key+' = '+str(params[key])+'\n')
    f.write('script run time = '+time.strftime('%H:%M:%S', time.gmtime(time.time()-start))+'\n\n')