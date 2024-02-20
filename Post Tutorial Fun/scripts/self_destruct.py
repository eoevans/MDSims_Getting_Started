import subprocess
import sys

if __name__ == '__main__':
    subprocess.run(['qdel',sys.argv[1]])