#import os
from os import environ

def runsdir():
    # RUNS = directory where run directories are stored
    # The hierarchy is such that $RUNS/$CODE/$EQUIL_NAME'_'$NUM contains the input and output of the run
    # Example : $RUNS/stella/w7xr003_0125
    return environ['RUNS']

def equilsdir():
    # EQUIL = directory where VMEC equilibria are stored
    # The hierarchy is sucha that $EQUIL/$EQUIL_NAME has inside the vmec .nc output and, when necessary
    # for Euterpe, also the vmec.equ.xxxx mapping files.
    # Example : $EQUIL/w7xr003 corresponds to the W7-X equilibrium ref. 3
    return environ['EQUIL']

