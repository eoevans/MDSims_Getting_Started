import MDAnalysis as mda
from MDAnalysis.transformations import translate
import sys
from params_validate import *
import numpy as np
from params_validate import get_surface_dims

def shift_box(input,output,dims):
    # shift surface coords so that the min coord is (0, 0, 0)
    # and the max coord is (x_len - dx, y_len - dy, z_len - dz)
    u = mda.Universe(input,transformations=translate(dims))
    u.atoms.write(output)

def shift_by(dims):
    return np.array(dims)/2

if __name__ == '__main__':
    # sys.argv[1] is the path to the pdb output by CHARMM-GUI
    # sys.argv[2] is the path to the pdb output by this script
    # sys.argv[3] is the path to the str output by CHARMM-GUI
    shift_box(sys.argv[1], sys.argv[2], shift_by(get_surface_dims(sys.argv[3])))