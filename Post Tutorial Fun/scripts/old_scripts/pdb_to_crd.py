import MDAnalysis as mda
from params_validate import read_config
import sys

def pdb_to_crd(input, output):
    u = mda.Universe(input)
    crd = mda.coordinates.CRD.CRDWriter(output)
    crd.write(u.atoms)

if __name__ == '__main__':
    files = read_config(sys.argv[1], 'files')
    pdb_to_crd(files['pdb_box'], files['crd'])