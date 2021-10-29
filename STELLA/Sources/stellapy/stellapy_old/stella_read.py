import numpy as np
from stella_dirs import *
from scipy.io import netcdf
#plt.rcParams.update({'font.size': 28})
#plt.rcParams['lines.linewidth'] = 2
import tabCompleter
from tabCompleter import *
from plotbox import *
from aux_functions import *
from os import listdir
from netCDF4 import *
import glob 
import os.path

# ==============================================================
# Some utils
def format1(value):
    return "%.3e" % value
def format2(value):
    return "%14.6e" % value
def format3(value):
    return "%4.2f" % value
def format4(value):
    return "%6.2f" % value
def format6(value):
    return "%7.3f" % value
def format5(value):
    return "%.5e" % value
def format7(value):
    return "%22.3f" % value
def format8(value):
    return "%04d" % value
def format9(value):
    return "%7.5f" % value
# Some utils ended
#===============================================================
   
def casestr(case=None):
    # Function that returns the string of the input, which
    # determines the name of the rest of output files.

    if case.endswith(".in"):
        buff  = case.split("/")
        return buff[size(buff)-1].split(".in")[0]
    else:
        if size(inputlist(case)) > 1:
            print("\nSpecify the input in the case field, more than one input file found:\n")
            print(inputlist(case))
            exit
        elif size(inputlist(case) == 1):
            return inputlist(case)[0].split(".in")[0]

def inputlist_r(case):
    inputs_level_0 = glob.glob(outdir(case)+'/*.in', recursive = True)
    inputs_level_1 = glob.glob(outdir(case)+'/*/*.in', recursive = True)
    return (inputs_level_0+inputs_level_1)
        
def inputlist(case, recursive=False):
    # Function that returns all the input file names
    # with extention ".in"
    inlist = []
    
    if recursive:
        inlist = inputlist_r(case=case)
    else:
        for f in listdir(outdir(case)):
            if f.endswith('.in'):
                if not f.startswith('.'):
                    inputname=f
                    inlist.append(f)
            
    return inlist

def outdir(case=None):
    if case.endswith(".in"):
        vcase=case.split("/")
        return runsdir()+'/'+ case.replace("/"+vcase[size(vcase)-1], '')
    else:
        return runsdir()+'/'+ case

def geotxtfile(case=None):
    # It returns the full path of an output file, endind with
    # the string value of "quant".
    if os.path.isfile(case):
        return case.split('.in')[0] + '.geometry'
    else:
        return outdir(case) + '/' + casestr(case) + '.geometry'

    
def outfile(case=None, quant=None):
    # It returns the full path of an output file, endind with
    # the string value of "quant".
    if os.path.isfile(case):
        return case.split('.in')[0] + '.' + quant
    else:
        return outdir(case) + '/' + casestr(case) + '.' + quant

def infile(case=None):
    # infile = input("Path to netcdf file: ")
    return outfile(case, quant='out.nc')

def fluxes_txt(case=None):
    # infile = input("Path to netcdf file: ")
    return outfile(case, quant='fluxes')

# ==================================================================
# Reading variables in the input *.in file

def torflux(case):
    # get torflux from input file.
    myfile  = open(outfile(case, quant='in'))
    content = float(myfile.read().split('torflux')[1].split('\n')[0].split('=')[1])
    return content

# ==================================================================
# Translation of quantities in stella_data module by Michael into
# functions with the run directory ("case") as single argument.
def read_stella_float(case, var):

    import numpy as np
   
    ncfile = netcdf.netcdf_file(infile(case),'r')
    
    try:
        arr = np.copy(ncfile.variables[var][:])
        flag = True
    except KeyError:
        print('INFO: '+var+' not found in netcdf file')
        arr = np.arange(1,dtype=float)
        flag = False
        
    return arr


def read_stella_value(case, var):
    woutfile = infile(case)
    d        = Dataset(woutfile, mode='r')
    return d.variables[var][:]

def kx(case):
    # get kx grid
    # this is the index of the first negative value of kx
    # note stella orders kx as (0, dkx, ..., kx_max, -kx_max, -kx_max+dkx, ..., -dkx)
    ncfile    = netcdf.netcdf_file(infile(case),'r')
    kx_stella = np.copy(ncfile.variables['kx'][:])
    nakx      = ncfile.dimensions['kx']
    nakx_mid  = nakx//2+1
    kx        = np.concatenate((kx_stella[nakx_mid:],kx_stella[:nakx_mid]))
    return kx, nakx, nakx_mid

def kx_stella(case):
    ncfile    = netcdf.netcdf_file(infile(case),'r')
    kx_stella = np.copy(ncfile.variables['kx'][:])
    return kx_stella

def ky(case):
    # get ky grid
    ncfile    = netcdf.netcdf_file(infile(case),'r')
    ky        = np.copy(ncfile.variables['ky'][:])
    naky      = ncfile.dimensions['ky']
    return ky, naky

def zed(case):
    # get zed grid
    ncfile    = netcdf.netcdf_file(infile(case),'r')
    zed       = np.copy(ncfile.variables['zed'][:])
    nzed      = zed.size
    iz0       = nzed//2+1
    return zed, nzed, iz0

def time(case):
    # get time grid
    ncfile    = netcdf.netcdf_file(infile(case),'r')
    time      = np.copy(ncfile.variables['t'][:])
    ntime     = time.size
    return time, ntime

def nspec(case):
    # number of kinetic species
    ncfile    = netcdf.netcdf_file(infile(case),'r')
    nspec     = ncfile.dimensions['species']
    return nspec

def geo(case):
    # get geometric quantities
    d         = Dataset(infile(case), mode='r')
    ncfile    = netcdf.netcdf_file(infile(case),'r')
    bmag      = np.copy(ncfile.variables['bmag'][:])
    gradpar   = np.copy(ncfile.variables['gradpar'][:])
    gbdrift   = np.copy(ncfile.variables['gbdrift'][:])
    gbdrift0  = np.copy(ncfile.variables['gbdrift0'][:])
    cvdrift   = np.copy(ncfile.variables['cvdrift'][:])
    cvdrift0  = np.copy(ncfile.variables['cvdrift0'][:])
    gds2      = np.copy(ncfile.variables['gds2'][:])
    gds21     = np.copy(ncfile.variables['gds21'][:])
    gds22     = np.copy(ncfile.variables['gds22'][:])
    shat      = float(d.variables['shat'][:])
    
    return bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, gds2, gds21, gds22, shat

def phi2_vs_kxky(case):
    # electrostatic potential averaged over z as function of (ky,kx,t)
    phi2_vs_kxky_stella = read_stella_float(case, 'phi2_vs_kxky')

#    phi2_vs_kxky_stella[:,0,0] = 0.0
#    phi2_vs_kxky = np.concatenate((phi2_vs_kxky_stella[:, kx(case)[2]:,:],\
#                                   phi2_vs_kxky_stella[:,:kx(case)[2] ,:]),axis=1)
    return phi2_vs_kxky_stella

def pflux_vs_kxky(case):
    pflux_vs_kxky_stella = read_stella_float(case, 'pflx_kxky')
    
    return pflux_vs_kxky_stella
    
def vflux_vs_kxky(case):
    vflux_vs_kxky_stella = read_stella_float(case, 'vflx_kxky')
    
    return vflux_vs_kxky_stella

def qflux_vs_kxky(case):
    qflux_vs_kxky_stella = read_stella_float(case, 'qflx_kxky')
    
    return qflux_vs_kxky_stella

def density_vs_kxky(case):
    density_vs_kxky_stella = read_stella_float(case, 'density')
    return density_vs_kxky_stella

def upar_vs_kxky(case):
    upar_vs_kxky_stella = read_stella_float(case, 'upar')
    return upar_vs_kxky_stella

def temperature_vs_kxky(case):
    temperature_vs_kxky_stella = read_stella_float(case, 'temperature')
    return temperature_vs_kxky_stella

def phi_vs_t(case):
    # electrostatic potential as a function of (z,kx,ky,t)
    phi_vs_t_stella = read_stella_float(case, 'phi_vs_t')
    return phi_vs_t_stella

def gvmus(case):
    # |g|^2 averaged over kx, ky, and z
    return  read_stella_float(case, 'gvmus')

def gzvs(case):
    # |g|^2 averaged over kx, ky, and mu
    return read_stella_float(case, 'gzvs')
    
def jacob(case):
    # jacobian for transformation to (rho,alpha,z) coordinates
    return read_stella_float(case, 'jacob')

def jtwist(case):
    # jtwist factor for twist-and-shift BC
    return read_stella_value(case, 'jtwist')

def grho(case):
    # gradient of normalized radial coordinate rho
    return read_stella_float(case, 'grho')

def phi2_stella(case):
    # modulus squared of electrostatic potential (averaged over space)
    return read_stella_float(case, 'phi2')

def es_part_flux(case):
    # time-dependent electrostatic particle flux for each species
    return read_stella_float(case, 'es_part_flux')

def es_heat_flux(case):
    # electrostatic heat flux
    return read_stella_float(case, 'es_heat_flux')

def es_mom_flux(case):
    # electrostatic momentum flux
    return read_stella_float(case, 'es_mom_flux')

def es_energy_exchange(case):
    return read_stella_float(case, 'es_energy_exchange')

def es_part_by_k(case):
    # time-dependent particle flux for each species as a function of (kx,ky)

    es_part_by_k_stella, es_part_by_k_present = \
                         read_stella_float(case, 'es_part_by_k')

    if es_part_by_k_present is not True:
        es_part_by_k_stella, es_part_by_k_present = \
                             read_stella_float(case, 'es_part_flux_by_mode')
        
    return es_part_by_k_stella, es_part_by_k_present

def es_mom_by_k(case):
    # time-dependent momentum flux for each species as a function of (kx,ky)
    es_mom_by_k_stella, es_mom_by_k_present = \
                        read_stella_float(case, 'es_mom_by_k')
    if es_mom_by_k_present is not True:
        es_mom_by_k_stella, es_mom_by_k_present = \
                            read_stella_float(case, 'es_mom_flux_by_mode')
    return es_mom_by_k_stella, es_mom_by_k_present

def es_energy_exchange_by_k(case):
    es_energy_exchange_by_k_stella, es_energy_exchange_by_k_present = \
                                    read_stella_float(case, 'es_energy_exchange_by_k')
    if es_energy_exchange_by_k_present is not True:
        es_energy_exchange_by_k_stella, es_energy_exchange_by_k_present = \
                                        read_stella_float(case, 'es_energy_exchange_by_mode')
    return es_energy_exchange_by_k_stella, es_energy_exchange_by_k_present

def es_energy_exchange_by_ky(case):
    return read_stella_float(case, 'es_energy_exchange_by_ky')

def vpa(case):
    # parallel velocity grid
    return read_stella_float(case, 'vpa')

def mu(case):
    # mu grid
    return read_stella_float(case, 'mu')

def es_part_sym(case):
    # electrostatic particle flux as function of (vpa,z)
    return read_stella_float(case, 'es_part_sym')

def es_heat_sym(case):
    # electrostatic heat flux as function of (vpa,z)
    return read_stella_float(case, 'es_heat_sym')

def es_mom_sym(case):
    # electrostatic momentum flux as function of (vpa,z)
    es_mom_sym_stella, es_mom_sym_present = read_stella_float(case, 'es_mom_sym')
    if vpa(case)[1] == False:
        es_mom_sym_present = False
    return es_mom_sym_stella, es_mom_sym_present

def xgrid(case):
    xgrid_stella, xgrid_present = \
           read_stella_float(case, 'xgrid')
    xgrid = np.concatenate((xgrid_stella[kx_stella(case).shape[0]//2+1:],\
                            xgrid_stella[:kx_stella(case).shape[0]//2+1]))
    return xgrid, xgrid_present

def dens(case):
    dens=read_stella_float(case, 'dens')
    dens_exp=factormult(dens,1e19)
    return dens_exp, size(dens)

def upar(case):
    # parallel flow fluctuation (kx,ky,z,t)
    return read_stella_float(case,'upar')

def temp(case):
    # temperature fluctuation (kx,ky,z,t)
    temp=read_stella_float(case,'temp')
    temp_exp=factormult(temp,1000)
    return temp_exp, size(temp)

def species(case):
    species=read_stella_float(case,'type_of_species')
    return species, size(species)

def nprim(case):
    return read_stella_float(case,'fprim')

def tprim(case):
    return read_stella_float(case,'tprim')

def charge(case):
    charge=read_stella_float(case,'charge')
    return charge, size(charge)

def mass(case):
    charge=read_stella_float(case,'mass')
    return charge, size(mass)

# ==================================================================
