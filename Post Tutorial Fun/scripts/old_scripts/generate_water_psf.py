import sys

def replace_in_str(str, replacement, index, start_at=True):
    if start_at:
        return str[:index]+replacement+str[index+len(replacement):]
    else:
        return str[:index-len(replacement)]+replacement+str[index:]

def generate(file, num_mols):
    OH2 = '       1 A    1    TP4E OH2  OT     0.000000       15.9994           0\n'
    OM = '       2 A    1    TP4E OM   LP    -1.048440        0.0000           0\n'
    H1 = '       3 A    1    TP4E H1   HT     0.524220        1.0080           0\n'
    H2 = '       4 A    1    TP4E H2   HT     0.524220        1.0080           0\n'
    bond = '       1       3       1       4       1       2       3       4\n'
    angle = '       3       1       2       3       1       4       4       1       2\n'
    with open(file, 'w') as f:
        f.write('PSF\n         3 !NTITLE\n')
        f.write('    '+str(num_mols*4)+' !NATOM\n')
        for i in range(num_mols):
            OH2 = replace_in_str(OH2, str(i*4+1), 8, start_at=False)
            OM = replace_in_str(OM, str(i*4+2), 8, start_at=False)
            H1 = replace_in_str(H1, str(i*4+3), 8, start_at=False)
            H2 = replace_in_str(H2, str(i*4+4), 8, start_at=False)
            OH2 = replace_in_str(OH2, str(i+1), 14, start_at=True)
            OM = replace_in_str(OM, str(i+1), 14, start_at=True)
            H1 = replace_in_str(H1, str(i+1), 14, start_at=True)
            H2 = replace_in_str(H2, str(i+1), 14, start_at=True)
            f.write(OH2)
            f.write(OM)
            f.write(H1)
            f.write(H2)
        f.write('\n    '+str(num_mols*4)+' !NBOND: bonds\n')
        for i in range(num_mols):
            bond = replace_in_str(bond, str(i*4+1), 8, start_at=False)
            bond = replace_in_str(bond, str(i*4+3), 16, start_at=False)
            bond = replace_in_str(bond, str(i*4+1), 24, start_at=False)
            bond = replace_in_str(bond, str(i*4+4), 32, start_at=False)
            bond = replace_in_str(bond, str(i*4+1), 40, start_at=False)
            bond = replace_in_str(bond, str(i*4+2), 48, start_at=False)
            bond = replace_in_str(bond, str(i*4+3), 56, start_at=False)
            bond = replace_in_str(bond, str(i*4+4), 64, start_at=False)
            f.write(bond)
        f.write('\n    '+str(num_mols)+' !NTHETA: angles\n')
        for i in range(num_mols//3):
            angle = replace_in_str(angle, str(i*12+3), 8, start_at=False)
            angle = replace_in_str(angle, str(i*12+1), 16, start_at=False)
            angle = replace_in_str(angle, str(i*12+4), 24, start_at=False)
            angle = replace_in_str(angle, str(i*12+7), 32, start_at=False)
            angle = replace_in_str(angle, str(i*12+5), 40, start_at=False)
            angle = replace_in_str(angle, str(i*12+8), 48, start_at=False)
            angle = replace_in_str(angle, str(i*12+11), 56, start_at=False)
            angle = replace_in_str(angle, str(i*12+9), 64, start_at=False)
            angle = replace_in_str(angle, str(i*12+12), 72, start_at=False)
            f.write(angle)
        f.write('\n       0 !NPHI: dihedrals\n\n       0 !NIMPHI: impropers\n\n       0 !NDON: donors\n\n       0 !NACC: acceptors\n\n')
        f.write('       0 !NNB\n\n       1       0 !NGRP\n       0       0       0')

if __name__ == '__main__':
    generate(sys.argv[1], int(sys.argv[2]))