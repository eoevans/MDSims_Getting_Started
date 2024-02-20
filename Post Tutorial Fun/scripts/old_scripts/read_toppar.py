from openmm import *
from openmm.app import *
from simtk.unit import *
import os

def read_toppar(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    for line in lines:
        