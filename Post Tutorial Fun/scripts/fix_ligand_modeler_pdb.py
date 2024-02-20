import sys

def fix_pdb(input_pdb, output_pdb, resname):
    with open(input_pdb, 'r') as f:
        lines = f.readlines()
    for i in range(len(lines)):
        if lines[i][:4] == 'ATOM':
            lines[i] = 'HETATM'+lines[i][6:17]+resname+lines[i][17+len(resname):]
    with open(output_pdb, 'w') as f:
        for line in lines:
            f.write(line)

if __name__ == '__main__':
    fix_pdb(sys.argv[1], sys.argv[2], sys.argv[3])