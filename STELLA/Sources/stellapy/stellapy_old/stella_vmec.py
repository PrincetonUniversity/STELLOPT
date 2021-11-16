import numpy
import scipy
import os

from os import *
from numpy import *
from scipy import interpolate
from netCDF4 import Dataset
from plotbox import *



def read_vmec_dataset(equil=None):   
    dataset = Dataset(woutfile(equil), mode='r')
    return dataset

def read_vmec_var(equil=None, varname='Aminor_p'):
    dataset = read_vmec_dataset(equil)
    return dataset.variables[varname][:]

def woutfile(equil=None, ext='nc'):
    equildir = equilsdir() + '/' + equil
    for f in listdir(equildir):
        if f.endswith("." + ext):
            woutfile = f
    return equildir + '/' + woutfile


def iota_stella(equil=None, svalue=None, plot=None):
    
    iotaf = read_vmec_var(equil, 'iotaf')
    iotas = read_vmec_var(equil, 'iotaf')
    ns    = read_vmec_var(equil, 'ns')
    svec  = linspace(0,1,ns)
    if svalue != None:
        sval    = svalue
        iotaval = interpol(svec, iotas, svalue)
    else:
        sval    = svec
        iotaval = iotas
        
    return sval, iotaval


# Some utilities and definitions necessary to avoid dependencies on
# other libraries.

def interpol(vector_x, vector_y, value_x, der=0):
    tck     = interpolate.splrep(vector_x, vector_y, s=0)
    value_y = interpolate.splev(value_x,tck,der=der)
    return value_y

def RZsym(R0=1.0, r=0.3, ntheta=100, delta=0.0, kappa=1.0, plot=False):
    theta = linspace(-pi, pi, ntheta)
    R     = R0+r*cos(theta+sin(theta)*arcsin(delta))
    Z     = kappa*r*sin(theta)
    
    if plot == True:
        ax    = pl2d(xrange=[0,R0+r], yrange=[-r,r], xlabel='$R$', ylabel='Z',\
                     fig_size=(8.5, 7.5), title=None, ax=None)
        ax.plot(R,Z,"-", color="black", linewidth=2)
        ax.plot(R0,0,"x", color="black", linewidth=2, markersize=10)
        ax.axis('equal')
        plt.show()
        
    return R, Z
    
def runsdir():
    # RUNS = directory where run directories are stored
    # The hierarchy is such that $RUNS/$CODE/$EQUIL_NAME'_'$NUM contains the input and output of the run
    # Example : $RUNS/stella/w7xr003_0125
    return os.environ['RUNS']

def equilsdir():
    # EQUIL = directory where VMEC equilibria are stored
    # The hierarchy is sucha that $EQUIL/$EQUIL_NAME has inside the vmec .nc output and, when necessary
    # for Euterpe, also the vmec.equ.xxxx mapping files.
    # Example : $EQUIL/w7xr003 corresponds to the W7-X equilibrium ref. 3
    return os.environ['EQUIL']


