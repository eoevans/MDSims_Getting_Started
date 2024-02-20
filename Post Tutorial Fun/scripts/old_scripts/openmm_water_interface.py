# This script was generated by OpenMM-Setup on 2023-02-28.

import sys
from openmm import *
from openmm.app import *
from simtk.unit import *
import time
import configparser
from params_validate import *
import os
from openmmtools import testsystems

def get_cutoff_distance(dims):
    # args are box dimensions in A. Returns half the box distance in nm or 1.4 nm per Singer
    cutoff_distance = 1.4
    for i in dims:
        if i/20 < cutoff_distance:
            cutoff_distance = i/20
    return cutoff_distance*nanometers

if __name__ == '__main__':
    # start the clock
    start = time.time()

    # read in parameters
    files, params = get_files_params(sys.argv[1], 'openmm')

    # get correct surface_dims
    surface_dims = get_surface_dims(files['surface_details'])
    # box_dims is in nanometers
    box_dims = [params['x_len'], params['y_len'], surface_dims[2]+params['z_len']]

    water_box = testsystems.FourSiteWaterBox(box_edge=3.0*nanometers, cutoff=1.0*nanometers, switch_width=None)
    forcefield = ForceField('charmm36.xml', 'charmm36/tip4pew.xml')
    psf = CharmmPsfFile(files['psf_surface'])
    pdb = PDBFile(files['pdb_surface'])

    modeller = Modeller(water_box.topology, water_box.positions)
    modeller.add(psf.topology, pdb.positions)    
    #ff = CharmmParameterSet(files['str'])

    # System Configuration

    nonbondedMethod = PME
    #nonbondedCutoff = 1.4*nanometers
    #nonbondedCutoff = 0.7*nanometers
    #ewaldErrorTolerance = 0.0005
    constraints = None
    rigidWater = False

    # Integration Options

    dt = 0.001*picoseconds
    friction = 1.0/picosecond
    pressure = 1.0*atmospheres
    barostatInterval = 25

    # Simulation Options
    #params['temp_step_size'] = (params['initial_temp'] - params['final_temp'])/params['temp_steps']
    #params['steps_per_temp'] = int(params['equil_steps']/params['temp_steps'])

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
    dataReporter = StateDataReporter(files['openmm_log'], 1000, totalSteps=params['equil_steps'],
        step=True, time=True, speed=True, progress=True, elapsedTime=True, remainingTime=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, separator='\t')

    # Prepare the Simulation

    print('Building system...')
    modeller.topology.setUnitCellDimensions([(params['x_len']/10), (params['y_len']/10), ((params['z_len']+surface_dims[2])/10)])
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=nonbondedMethod,
                              nonbondedCutoff=get_cutoff_distance([params['x_len'],params['y_len'],params['z_len']]), 
                              constraints=constraints, rigidWater=rigidWater) 

    # system.addForce(MonteCarloBarostat(pressure, temperature, barostatInterval))
    #integrator = LangevinMiddleIntegrator(params['initial_temp']*kelvin, friction, dt)
    integrator = VerletIntegrator(dt)
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    #simulation.context.setPeriodicBoxVectors((params['x_len']/10)*nanometers, (params['y_len']/10)*nanometers, ((params['z_len']+surface_dims[2])/10)*nanometers)
    # print(f"inpcrd.boxVectors is {inpcrd.boxVectors}\n")
    #if inpcrd.boxVectors is not None:
    #    print(f"inpcrd.boxVectors is {inpcrd.boxVectors}\n")
    #simulation.context.setPeriodicBoxVectors(Vec3(params['side_len']/10*nanometers, params['side_len']/10*nanometers, params['side_len']/10*nanometers))

    # Minimize and Equilibrate

    print('Performing energy minimization...')
    simulation.minimizeEnergy()
    print('Equilibrating...')
    simulation.reporters.append(dataReporter)
    for i in range(int(params['equil_steps']/params['v_to_t_steps'])):
        simulation.step(params['v_to_t_steps'])
        simulation.context.setVelocitiesToTemperature(params['temp'])


    # Simulate
    print('Simulating...')
    simulation.reporters.append(dcdReporter)
    simulation.reporters.append(PDBReporter('MDSim/data/simulation_results/water/tip4p2005/run1/output.pdb', 1000))
    # simulation.reporters.append(pdbReporter)
    simulation.currentStep = 0
    simulation.step(params['sim_steps'])


    # end simulation section -- now to write the simulation parameters to the results file
    params['num_trajectory_frames'] = params['sim_steps']/dcd_interval
    params['run_time'] = params['sim_steps']*dt
    params['frame_spacing'] = params['run_time']/params['num_trajectory_frames']
    params['equil_time'] = params['equil_steps']*dt

    with open(files['results'], 'w') as f:
        f.write('Openmm\n')
        f.write('All temperatures are in Kelvin.\n')
        for key in params:
            f.write(key+' = '+str(params[key])+'\n')
        f.write('script run time = '+time.strftime('%H:%M:%S', time.gmtime(time.time()-start))+'\n\n')