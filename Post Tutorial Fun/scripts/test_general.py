from CalcT12_anisotropic import values_of_interest
import numpy as np
import sys
from params_validate import *
import subprocess
import time

if __name__ == '__main__':
    submit_fast_job('MDSim/scripts/self_destruct.py', 'MDSim/submit/self_destruct.sh', get_job_number('test'))
    time.sleep(600)