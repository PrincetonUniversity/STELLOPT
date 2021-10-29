import os
import matplotlib.pyplot as plt 

def load_styleFigures():
    plt.rcParams['lines.linewidth'] = 2 
    plt.rc('font', family='serif')
    plt.rc('font', size=20)
    plt.rcParams.update({'figure.autolayout': True})
    
    # Marconi does not support the fancy latex font
    divider = '\\' if (os.name == 'nt') else '/'
    if os.getcwd().split(divider)[1] == 'marconi_work':
        plt.rc('text', usetex=False)
    if os.getcwd().split(divider)[1] == 'home':
        plt.rc('text', usetex=True)

# Force it to load before opening a pyplot interactive loop.
load_styleFigures()
