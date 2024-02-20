import sys
import os
from MDSim_automate import get_mols
from params_validate import read_config

def fix_gaussian_pdb(pdb_old, pdb_middle, pdb_new, name):
    with open(pdb_old, 'r') as f:
        lines = f.readlines()
    with open(pdb_middle, 'w') as f:
        f.write(lines[0][:11]+name+'\n')
        f.write(lines[1])

        i = 2
        while lines[i] != "END\n":
            f.write(lines[i][:17]+name+lines[i][17+len(name):])
            i += 1
        while i < len(lines):
            f.write(lines[i])
            i += 1
        f.close()
    os.system('pdb_sort '+pdb_middle+' | pdb_element | pdb_tidy | pdb_b -10 | pdb_occ -1 | pdb_uniqname > '+pdb_new)

if __name__ == '__main__':
    # add the residue name column to a gaussian generated pdb, and prepare it for use with packmol using pdb_tools
    # first argument should be the gaussian pdb path, titled molecule.pdb
    # second argument should be the named pdb path, titled molecule_named.pdb
    # the third argument should be the valid pdb path, titled molecule_valid.pdb
    # the fourth argument is the residue name. This should match the residue name in the str/rtf file
    fix_gaussian_pdb(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])