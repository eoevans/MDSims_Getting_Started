import matplotlib.pyplot as plt
import os
import sys

if __name__ == '__main__':
    files = os.listdir(sys.argv[1])
    params = [20, 50, 100, 200, 400, 800]
    vals = [[], [], []]
    count = 0
    for i in files:
        
        if i.split('.')[1][-5:] == 'indiv':
            with open (i, 'r') as f:
                lines = f.readlines()
            for i in range(2, len(lines)):
                if (i-2) % 4 == 0:
                    pass
                else:
                    line = lines[i].split()
                    vals[]

