import MDAnalysis as mda
import sys
from params_validate import get_files_params
from make_box import get_pdbs

def pdb_merge(pdb_surface, pdbs_liquid, merged_pdb):
    # combine many pdbs into one.
    # used for making the mixed surface/liquid pdb. I tried doing it with vmd in the same step as generating the psf
    # but those pdbs were giving errors in openmm
    u_final = mda.Universe(pdb_surface)
    for i in pdbs_liquid:
        u_add = mda.Universe(i)
        u_final = mda.Merge(u_final.atoms, u_add.atoms)
    u_final.atoms.write(merged_pdb)

if __name__ == '__main__':
    files, params = get_files_params(sys.argv[1], 'packmol_vmd')
    pdb_merge(files['pdb_surface'], get_pdbs(params, files, inp=False), files['pdb_sim'])
    #pdb_merge(sys.argv[1], [sys.argv[2]], sys.argv[3])