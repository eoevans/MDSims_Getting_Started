import math
from params_validate import validate_section

pi = math.pi
kb = 1.380649*(10**-23)
squiggle = 2.837297
Na = 6.02214*(10**23)

# when simulating a new molecule for the first time, the references in shear_vis, count_mols_in_box, get_name, get_CHARMM_files must be updated
def shear_vis(mol):
    # ref_mols here contains the shear viscosity of each molecule. This has units of cP or mPa*s
    ref_mols = {'C5H12': 0.24,'C6H14': 0.31,'C7H16': 0.42,'C8H18': 0.55,'C10H22': 0.93,'C12H26': 1.51,
                'C14H30': 2.34,'C16H34': 3.48,'C18H38': None,'C20H42': None,'acetone': 0.32,'ethanol': 1.074,
                'glycerol': 954, 'isopropanol': 2.038, 'water': 1}
    for key in ref_mols:
        if mol == key:
            return ref_mols[key]/1000

def count_mols_in_box(mol, dims):
    # density is in g/mL
    ref_density = {'C5H12': 0.626,'C6H14': 0.659,'C7H16': 0.684,'C8H18': 0.703,'C10H22': 0.730,'C12H26': 0.750,
                    'C14H30': 0.762,'C16H34': 0.773,'C18H38': 0.777,'C20H42': 0.789,'acetone': 0.791,'ethanol': 0.79,
                    'toluene': 0.865, '1-2-diaminopropane': 0.8584, 'DGEBA': 1.16, 'ethylene_glycol': 1.115, 'glycerol': 1.261,
                    'isopropanol': 0.785, 'water': 1}
    # molar mass is in g/mol
    ref_molar_mass = {'C5H12': 72.15,'C6H14': 86.18,'C7H16': 100.21,'C8H18': 114.23,'C10H22': 142.29,'C12H26': 170.33,
                    'C14H30': 198.39,'C16H34': 226.41,'C18H38': 254.494,'C20H42': 282.5475,'acetone': 58.08, 'ethanol': 46.07,
                    'toluene': 92.14, '1-2-diaminopropane': 74.13, 'DGEBA': 340.4, 'ethylene_glycol': 62.07, 'glycerol': 92.09,
                    'isopropanol': 60.10, 'water': 18.02}
    for key in ref_density:
        if mol == key:
            return round((ref_density[key]/ref_molar_mass[key]*Na)*((dims[0]*dims[1]*dims[2])*10**-24))
        
def get_segment_name(file):
    ref_dirs = {'MDSim/data/surface_data/aluminum': 'ALM', 'MDSim/data/surface_data/alumina': 'AO00', 'MDSim/data/surface_data/graphite': 'NM'}
    for key in ref_dirs:
        if file[:len(key)] == key:
            return ref_dirs[key]
        
def get_single_run_time(section_str, sim_time, equil_time, box_volume):
    # sent the section for which the time is calculated, the equilibration simulation time of the openmm run
    # and the box volume including the surface volume, return the approximate run time
    # this run time is scaled by the number of mols and runs in MDSim_automate
    ref_constants = {'packmol_vmd_sub': 20, 'openmm_sub': 0.19, 'inter_sub': 0.0055,
                    'intra_sub': 0.18, 'results_sub': 1, 'diffusion_sub': 10, 'pdb_merge_sub': 1}
    if section_str == 'openmm_sub':
        return ref_constants[section_str]*box_volume*(sim_time + equil_time/2)
    elif section_str == 'inter_sub':
        return ref_constants[section_str]*(box_volume**2)*(sim_time)
    elif section_str == 'intra_sub':
        return ref_constants[section_str]*(box_volume*sim_time)
    else:
        # currently, time scaling based on openmm_time and box size is only implemented for the above sections
        # should also be implemented for packmol, which will scale with box size (not sure how), and diffusion
        # which sould scale with sim_time and box_size^2 (like inter).  
        return ref_constants[section_str]

def get_CHARMM_files(mol):
    ref_mols = {'C5H12': 'cgenff','C6H14': 'cgenff','C7H16': 'lipid','C8H18': 'cgenff','C10H22': 'lipid','C12H26': 'cgenff',
                'C14H30': 'lipid','C16H34': 'lipid','C18H38': 'cgenff','C20H42': 'cgenff','acetone': 'cgenff','ethanol': 'cgenff', 
                'toluene': 'cgenff', '1-2-diaminopropane': 'cgenff', 'DGEBA': 'cgenff', 'ethylene_glycol': 'cgenff', 'glycerol': 'cgenff',
                'isopropanol': 'cgenff', 'aluminum': 'interface', 'water': 'water'}
    for key in ref_mols:
        if key == mol:
            if ref_mols[key] == 'cgenff':
                return '\nstr = MDSim/toppar/'+mol+'.str\nrtf = MDSim/toppar/top_all36_cgenff.rtf\nprm = MDSim/toppar/par_all36_cgenff.prm\n'
            elif ref_mols[key] == 'lipid':
                return '\nstr = MDSim/toppar/toppar_all36_lipid_model.str\nrtf = MDSim/toppar/top_all36_lipid.rtf\nprm = MDSim/toppar/par_all36_lipid.prm\n'
            elif ref_mols[key] == 'interface':
                return '\nstr = MDSim/toppar/empty.str\nrtf = MDSim/toppar/top_interface.rtf\nprm = MDSim/toppar/par_interface.prm\n'
            elif ref_mols[key] == 'water':
                return '\nstr = MDSim/toppar/toppar_water_ions.str\nrtf = MDSim/toppar/top_all36_cgenff.rtf\nprm = MDSim/toppar/par_all36_cgenff.prm\n'