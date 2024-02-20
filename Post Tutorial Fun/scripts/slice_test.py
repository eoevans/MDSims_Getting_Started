import MDAnalysis as mda
import sys
from params_validate import read_config

files = read_config(sys.argv[1], 'files')
u = mda.Universe(files['pdb_sim'],files['dcd'])
combined = u.atoms
for ts in u.trajectory:
    in_box = u.select_atoms('prop z >= 62.15')
    combined = combined & in_box
print(len(combined)/(len(u.select_atoms('resname OCN'))/2))