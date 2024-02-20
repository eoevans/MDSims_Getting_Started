import MDAnalysis
from params_validate import read_config
import sys
from params_validate import *


files = read_config(sys.argv[1], 'files')
u = MDAnalysis.Universe(files['pdb_sim'], files['dcd'])
with MDAnalysis.Writer(files['dcd'].split('.')[0]+'_'+sys.argv[2]+'.dcd', u.atoms.n_atoms) as W:
    for ts in u.trajectory[:int(sys.argv[3])]:
        W.write(u.atoms)