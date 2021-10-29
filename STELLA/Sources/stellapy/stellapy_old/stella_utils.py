import h5py
from stella_dirs import *
from stella_vmec import *
from numpy import *
from pylab import *
from struct import *
from scipy import *
from matplotlib import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import LinearLocator, FixedLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
from scipy.io import netcdf
import io
plt.rcParams.update({'font.size': 28})
plt.rcParams['lines.linewidth'] = 2

import tabCompleter

from read_wout import *
from plotbox import *
from tabCompleter import *
import struct
import physcon as pc
import time
from stella_read import *


##def inputlist(case):
##    i      = 1
##    inlist = []

##    for f in listdir(outdir(case)):
##        if f.endswith(".in"):
##            inputname=f
##            i = i +1
##            inlist.append(f)
##            print(os.path.join(outdir(case), f))

##    if i > 1: morethanone = True
    
##    return inlist

def interpol(vector_x, vector_y, value_x, der=0):

    tck     = interpolate.splrep(vector_x, vector_y, s=0)
    value_y = interpolate.splev(value_x,tck,der=der)

    return value_y

def nfield_periods(equil, svalue=None, dtheta=3, wrt=False):
    # dtheta (in units of pi) is the range the field
    # flux tube will extent along poloidally. I.e. 
    # A flux tube centered at (theta, zeta) = ( iota*zeta_center , zeta_center) will
    # cover the range iota*zeta_center +/- 3pi poloidally. 
    nfp = read_vmec_var(equil, varname='nfp')
    return nfp * dtheta / iota(equil, svalue)[1]


def ref_values(case=None, verbose=False, code='stella', tref=None, nref=None):
    if tref == None:
        tref = temp(case)[0][0]   # in eV
    else:
        # For simplicity we introduce externally T in keV
        tref = tref * 1000.
        
    if nref == None:
        nref = dens(case)[0][0]   # in eV
    else:
        # For simplicity we introduce externally T in units of 10^19 m^-3
        nref = nref * 1E19
        
    tref = tref * pc.e        # in J
    mref = mass(case)[0][0]   # in units of m_p
    mref = mref * pc.m_p      # in kg
    zref = charge(case)[0][0] # in units of e
    qe   = pc.e

    # We read the first line of the txt vmec_geo file with ref values
    # l_r=a and B_ref
    vmec_geo = geotxtfile(case)
    a = io.open(str(vmec_geo),"r")
    a.readline()
    ref_geo_values = a.readline().split("#")[1].strip().split(" ")
    ref_geo_values = [float(value) for value in ref_geo_values if value != '']
    rhoc, qinp, shat, rhotor, lref, Bref, dxdpsi, dydalpha = ref_geo_values

    vth_ref = sqrt(2*tref/mref)
    if code == 'gene': vth_ref = sqrt(tref/mref)
    ome_ref = zref*pc.e*Bref/mref
    rho_ref = vth_ref / ome_ref 

    gamma_ref = (rho_ref/lref)**2.0 * nref * vth_ref
    mom_ref   = (rho_ref/lref)**2.0 * nref * lref * mref * vth_ref**2.0
    heat_ref  = (rho_ref/lref)**2.0 * nref * tref * vth_ref

    if verbose:
        print('')
        print('l_ref := a =', lref, 'm.')
        print('B_ref      =', Bref, 'T.')
        print('T_ref      =', tref, 'J. =', tref/pc.e, 'eV.')
        print('n_ref      =', nref, 'm^{-3}.')
        print('m_ref      =', mref, 'kg.')
        print('Z_ref      =', zref, 'adimensional.')
        if code == 'stella':
            print('vth_ref   := sqrt(2*T_ref/m_ref) =', vth_ref, 'm/s.')
        elif code == 'gene':
            print('cs_ref    := sqrt(T_ref/m_ref) =', vth_ref, 'm/s.')
        print('omega_ref := eZ_ref*B_ref/m_ref  =', ome_ref, 's^{-1}.')
        print('rho_ref   := vth_ref/rho_ref     =', rho_ref, 'm.')
        if code == 'stella':
            print('Q_ref     := ((rho_ref/l_ref)**2) * n_ref * T_ref * vth_ref =', heat_ref, 'W/m^2')
        elif code == 'gene':
            print('Q_ref     := ((rho_ref/l_ref)**2) * n_ref * T_ref * cs_ref =', heat_ref, 'W/m^2')
        print('')
    
    return tref, nref, mref, zref, vth_ref, ome_ref, rho_ref,\
           lref, Bref, gamma_ref, mom_ref, heat_ref


def merge(case, quant='phi2_vs_kxky'):
    # This functions merge the different arrays of a certain
    # quantity and retunrs it.
    # It assumes that the input files that are present in
    # outdir(case) and all its subdirectories correspond
    # to a run with multiple restarts.
    input_files = inputlist(case=case, recursive=True)
    dim         = size(input_files)
    time_merged = []
    t0, tf      = empty(dim, dtype='float'), empty(dim, dtype='float')

    for i in input_files:
        t0[input_files.index(i)]=time(i)[0][0]
        tf[input_files.index(i)]=time(i)[0][size(time(i)[0])-1]
        
    # We sort the arrays by increasing value of initial simulated time
    inds               = t0.argsort()
    input_files        = list(array(input_files)[inds])
    t0                 = t0[inds]
    tf                 = tf[inds]

    time_merged = time(input_files[0])[0]
    if quant == 'fluxes':
        quantity_merged = loadtxt(fluxes_txt(input_files[0]))
    else:
        quantity_merged = read_stella_float(input_files[0], quant)
        

    for i in input_files[1:]:
            #time_merged = concatenate(time_merged, time(i)[0] > tf[input_files.index(i)-1])
            time_temp  = time(i)[0][time(i)[0] > tf[input_files.index(i)-1]]
            time_merged = np.concatenate((time_merged, time_temp), axis=0)
            
            if quant == 'fluxes':
                quantity_temp = loadtxt(fluxes_txt(i))[time(i)[0] > tf[input_files.index(i)-1],:]
                quantity_merged = np.concatenate((quantity_merged, quantity_temp), axis=0)
            else:
                quantity_temp = read_stella_float(i, quant)[time(i)[0] > tf[input_files.index(i)-1],:,:]
                quantity_merged = np.concatenate((quantity_merged, quantity_temp), axis=0)

    return(time_merged, quantity_merged)

def ngrid_to_nk(nx=20, ny=20):
    # This function answers:
    # How many modes (nkx,nky) are considered when the box
    # is devided in (nx,ny) segments?
    nkx = 2*((nx-1)/3) + 1
    nky = (ny-1)/3 + 1
    return(nkx, nky)

def nk_to_ngrid(nkx=20, nky=20):
    # This function answers:
    # How many modes (nkx,nky) are considered when the box
    # is devided in (nx,ny) segments?
    nx  = 3*((nkx - 1)/2) + 1
    ny  = 3*(nky - 1)+1
    return(nx, ny)


def exp_prof(file=None, fsl_col=0,  quant_col=1,
             fsl='rho', rho_val=0.5,smooth=0.0, ns=100, plot=False, a=None):

    data = loadtxt(file, dtype='float', usecols=(fsl_col, quant_col))

    # a/L_X = a * dX/dr = a * (drho/dr) * dlog(X)/drho = dlogX/drho
    # a/L_X = a * dX/dr = a * (ds/dr)   * dlog(X)/ds   = 2*sqrt(s) (dlogX/ds)

    fsl_vec = np.linspace(0,1.,ns)
    
    if fsl == 'rho':
        tck     = interpolate.splrep(data[:,0],data[:,1],s=smooth)
        # We build the profiles
        quant     = interpolate.splev(fsl_vec,tck,der=0)
        a_over_Lq = interpolate.splev(fsl_vec,tck,der=1)/quant
        
        # We get the values at the specified s (or psitor, or svalue)
        quant_at_loc     = interpolate.splev(rho_val,tck,der=0)
        a_over_Lq_at_loc = interpolate.splev(rho_val,tck,der=1)/quant_at_loc

        xlab = '$r/a$'
        
    if fsl == 's':
        tck     = interpolate.splrep(data[:,fsl_col],data[:,quant_col],s=smooth)
        # We build the profiles
        quant     = interpolate.splev(fsl_vec, tck,der=0)
        a_over_Lq = 2*sqrt(fsl_vec)*interpolate.splev(fsl_vec,tck,der=1)/quant
        
        # We get the values at the specified s (or psitor, or svalue)
        quant_at_loc     = interpolate.splev(rho_val**2.0,tck,der=0)
        a_over_Lq_at_loc = 2*rho_val*interpolate.splev(rho_val**2.0,tck,der=1)/quant_at_loc
        print(interpolate.splev(rho_val**2.0,tck,der=1)/quant_at_loc)

        xlab = '$s$'

    if fsl == 'r':
        data[:,fsl_col] = data[:,fsl_col]/a
        tck     = interpolate.splrep(data[:,fsl_col],data[:,quant_col],s=smooth)
        # We build the profiles
        quant     = interpolate.splev(fsl_vec, tck,der=0)
        a_over_Lq = interpolate.splev(fsl_vec, tck,der=1)/quant
        
        # We get the values at the specified s (or psitor, or svalue)
        quant_at_loc     = interpolate.splev(rho_val,tck,der=0)
        a_over_Lq_at_loc = interpolate.splev(rho_val,tck,der=1)/quant_at_loc

        xlab = '$r/a$'        
    # plotting

    if plot == True:
        ylab = 'Q'
        fig = plt.figure(figsize=(18, 9))
        fig.subplots_adjust(left=0.1, wspace=0.35)
        
        ax1 = fig.add_subplot(121)
        title = '$Q$ at $r/a = $' + str(rho_val)+ ') = ' + str(quant_at_loc)
        pl2d(xrange=[0,1], xlabel=xlab, ylabel=ylab,\
             fig_size=(8.5, 7.5), title=title, ax=ax1)
        ax1.autoscale()
        ax1.plot(data[:, 0], data[:, 1], 's', color='red',\
                 markeredgewidth=2, markersize=8.0, markerfacecolor='white', label='Exp.')
        ax1.plot(fsl_vec, quant, '-', color='navy',\
                 linewidth=4, markerfacecolor='white', label='Interpol.')
        ax1.legend(loc='best', labelspacing=0.1, prop={'size':20})
        
        title = '$a/L_Q$ at $r/a = $' + str(rho_val)+ '$ = $ ' + str(a_over_Lq_at_loc)
        ylab = '$a/L_Q$'
        ax2 = fig.add_subplot(122)
        pl2d(xrange=[0,1], xlabel=xlab, yrange=None, ylabel=ylab,\
             fig_size=(8.5, 7.5), title=title, ax=ax2)
        ax2.autoscale()
        ax2.plot(fsl_vec, a_over_Lq, '-', color='navy',\
                 linewidth=4, markerfacecolor='white', label='Interpol.')
        ax2.legend(loc='best', labelspacing=0.1, prop={'size':20})   
        show()

    return rho_val, quant_at_loc[()], -a_over_Lq_at_loc

def euprof(fileprof=None, svalue=None):
    # This function reads a profile in the format of EUTERPE
    # and converts it to another with the parameters
    # that stella needs.
    profdata=loadtxt(fileprof, dtype='float')
    if shape(profdata)[1] == 9:
        S,DTI,TI,DTE,TE,DNI,NI,DNE,NE = arange(0,9)

    dsdrho = 2*sqrt(profdata[:,S])

    s       =  profdata[:,S]
    nine    =  profdata[:,NI]/profdata[:,NE]
    tite    =  profdata[:,TI]/profdata[:,TE]
    tprim_i = -profdata[:,DTI]*dsdrho
    tprim_e = -profdata[:,DTE]*dsdrho
    fprim_i = -profdata[:,DNI]*dsdrho
    fprim_e = -profdata[:,DNE]*dsdrho
    t_i     =  profdata[:,TI]/1000.  # stella uses keV for temp
    t_e     =  profdata[:,TE]/1000.  # stella uses keV for temp
    dens_e  =  profdata[:,NE]/1.0E19 # stella uses keV for temp
    dens_i  =  profdata[:,NI]/1.0E19 # stella uses keV for temp

    if svalue != None:
        print("&vmec_parameters")
        print("torflux="+str(format9(svalue))+'\n')
        print("&parameters")
        print("nine="  +str(format9(interpol(profdata[:,S],nine,svalue))))
        print("tite="  +str(format9(interpol(profdata[:,S],tite,svalue)))+'\n')
        print("&species_parameters")
        print("dens="  +str(format9(interpol(profdata[:,S],dens_i,svalue))))
        print("temp="  +str(format9(interpol(profdata[:,S],t_i,   svalue))))
        print("tprim=" +str(format9(interpol(profdata[:,S],tprim_i[:],svalue))))
        print("fprim=" +str(format9(interpol(profdata[:,S],fprim_i[:],svalue)))+'\n')
        print("\n")
        print("&species_parameters_2")        
        print("dens="  +str(format9(interpol(profdata[:,S],dens_e,svalue))))
        print("temp="  +str(format9(interpol(profdata[:,S],t_e,   svalue))))
        print("tprim=" +str(format9(interpol(profdata[:,S],tprim_e,svalue))))
        print("fprim=" +str(format9(interpol(profdata[:,S],fprim_e,svalue))))

        
        

def h5_to_stella(equil=None, file=None, Ti_diag='CXRS',\
                  title=None, output='old', plot=True, svalue=None):

    # Ti: it can take XICS or CXRS
    #
    f           = h5py.File(file, 'r')
    filename    = file.split("/")[1].split(".h5")[0]

    print('\n')
    for i in list(f.keys()):
        print('Key = ', str(i), '\t ---> subkeys: ', list(f[i].keys()))
    print('\n')    

    if 'XICS_data' in f['Ti'].keys():
        xics = True
        xics = False
        Ti_exp_XICS = f['Ti']['XICS_data']
        Ti_fit_XICS = f['Ti']['XICS_fit']
    if 'CXRS_data' in f['Ti'].keys():
        cxrs = True
        Ti_exp_CXRS = f['Ti']['CXRS_data']
        Ti_fit_CXRS = f['Ti']['CXRS_fit']

    Te_exp_TS   = f['Te']['data']
    Te_fit_TS   = f['Te']['fit']
    ne_exp_TS   = f['ne20']['data']
    ne_fit_TS   = f['ne20']['fit']
    ni_fit_TS   = ne_fit_TS

    # Profile interpolating
    
    nine = ni_fit_TS[1,:]/ne_fit_TS[1,:]
    
    if Ti_diag == 'CXRS':
        tite    = Ti_fit_CXRS[1,:]/Te_fit_TS[1,:]
        t_i     = Ti_fit_CXRS[1,:]
        tprim_i = -Ti_fit_CXRS[2,:]/Ti_fit_CXRS[1,:]
    elif Ti_diag == 'XICS':
        tite = Ti_fit_XICS[1,:]/Te_fit_TS[1,:]
        t_i     = Ti_fit_XICS[1,:]
        tprim_i = -Ti_fit_XICS[2,:]/Ti_fit_CXRS[1,:]
        
    t_e     = Te_fit_TS[1,:]
    tprim_e = -Te_fit_TS[2,:]/Te_fit_TS[1,:]
    fprim_e = -ne_fit_TS[2,:]/ne_fit_TS[1,:]
    fprim_i = fprim_e
    rho     = Te_fit_TS[0,:]

    if svalue != None:
        print("&vmec_parameters")
        print("torflux="+str(format9(svalue))+'\n')
        print("&parameters")
        print("nine="  +str(format9(interpol(rho,nine,sqrt(svalue)))))
        print("tite="  +str(format9(interpol(rho,tite,sqrt(svalue))))+'\n')
        print("&species_parameters")
        # Remember stella nref is 1E19
        print("dens="  +str(format9(interpol(rho,ni_fit_TS[1,:] ,sqrt(svalue))*10)))
        print("temp="  +str(format9(interpol(rho,t_i     ,sqrt(svalue)))))
        print("tprim=" +str(format9(interpol(rho,tprim_i ,sqrt(svalue)))))
        print("fprim=" +str(format9(interpol(rho,fprim_i ,sqrt(svalue))))+'\n')
        print("&species_parameters_2")        
        print("dens="  +str(format9(interpol(rho,ne_fit_TS[1,:],sqrt(svalue))*10))) 
        print("temp="  +str(format9(interpol(rho,t_e    ,sqrt(svalue)))))
        print("tprim=" +str(format9(interpol(rho,tprim_e,sqrt(svalue)))))
        print("fprim=" +str(format9(interpol(rho,fprim_e,sqrt(svalue))))+'\n')


    # Profile plotting
    if plot == True:
        # Temperature
        ax = pl2d(xrange=[0,1], yrange=[0,3], xlabel='$r/a$', ylabel='$T_{e}, T_{i}$ [keV]',\
                  fig_size=(8.5, 7.5), title=title)
        if xics: 
            ax.errorbar(Ti_exp_XICS[0,:], Ti_exp_XICS[1,:], Ti_exp_XICS[2,:],\
                        linewidth=2.0, fmt='s', color='navy', markersize=10,\
                        mfc='white', markeredgewidth=2.0, mec='navy',\
                        label='$T_{i}^{\\mathrm{XICS}}$')
            ax.plot(Ti_fit_XICS[0,:], Ti_fit_XICS[1,:], '-', linewidth=3.0, color='navy',)

        if cxrs:
            ax.errorbar(Ti_exp_CXRS[0,:], Ti_exp_CXRS[1,:], Ti_exp_CXRS[2,:],\
                        linewidth=2.0, fmt='o', color='cornflowerblue', markersize=10,\
                        mfc='white', markeredgewidth=2.0, mec='cornflowerblue',\
                        label='$T_{i}^{\\mathrm{CXRS}}$')
            ax.plot(Ti_fit_CXRS[0,:], Ti_fit_CXRS[1,:], '-', linewidth=3.0, color='cornflowerblue')
        
        ax.errorbar(Te_exp_TS[0,:], Te_exp_TS[1,:], Te_exp_TS[2,:],\
                    linewidth=2.0, fmt='^', color='crimson', markersize=10,\
                    mfc='white', markeredgewidth=2.0, mec='crimson',\
                    label='$T_{e}^{\\mathrm{TS}}$')
        ax.plot(Te_fit_TS[0,:], Te_fit_TS[1,:], '-', linewidth=3.0, color='crimson')
        ax.legend(loc='best', labelspacing=0.1, prop={'size':20})
        
        # Density
        ax = pl2d(xrange=[0,1], yrange=[0,1], xlabel='$r/a$', ylabel='$n_e$ [$10^{20}$ m$^{-3}$]',\
                  fig_size=(8.5, 7.5), title=title)
        ax.errorbar(ne_exp_TS[0,:], ne_exp_TS[1,:], ne_exp_TS[2,:],\
                    linewidth=2.0, fmt='s', color='limegreen', markersize=10,\
                    mfc='white', markeredgewidth=2.0, mec='limegreen',\
                    label='$n_{i}^{\\mathrm{TS}}$')
        ax.plot(ne_fit_TS[0,:], ne_fit_TS[1,:], '-', linewidth=3.0, color='limegreen',)
        ax.legend(loc='best', labelspacing=0.1, prop={'size':20})
        
        # Gradient length-scales
        ax = pl2d(xrange=[0,1], yrange=[-10,3], xlabel='$r/a$', ylabel='$a/L_{X}$',\
                  fig_size=(8.5, 7.5), title=title)
        ax.plot(ne_fit_TS[0,:], ne_fit_TS[2,:]/ne_fit_TS[1,:], '-s', markersize=10,\
                mfc='white', markeredgewidth=2.0, mec='limegreen',\
                linewidth=3.0, color='limegreen',\
                label='$X=n_{e}^{\\mathrm{TS}}$')
        if cxrs:
            ax.plot(Ti_fit_CXRS[0,:], Ti_fit_CXRS[2,:]/Ti_fit_CXRS[1,:], '-o',\
                    color='cornflowerblue', markersize=10,\
                    mfc='white', markeredgewidth=2.0, linewidth=3.0,\
                    label='$X=T_{i}^{\\mathrm{CXRS}}$')
        if xics:
            ax.plot(Ti_fit_XICS[0,:], Ti_fit_XICS[2,:]/Ti_fit_XICS[1,:], '-s', linewidth=3.0, color='navy',\
                    markersize=10, mfc='white', markeredgewidth=2.0, mec='navy',\
                    label='$X=T_{i}^{\\mathrm{XICS}}$')
            
        ax.plot(Te_fit_TS[0,:], Te_fit_TS[2,:]/Te_fit_TS[1,:], '-^', markersize=10,\
                mfc='white', markeredgewidth=2.0, mec='crimson',\
                linewidth=3.0, color='crimson', label='$X=T_{e}^{\\mathrm{TS}}$')
        ax.legend(loc='best', labelspacing=0.1, prop={'size':20})
        plt.show()
        



    
##    f       = open('./stella_input_prof.dat',"w")
##    f.write('# (1)psi\t(2)nine\t\t(3)tite\t(4)dens(i)\t(5)temp(i)\t(6)tprim(i)\t(7)fprim(i)'+\
##            '\t(8)dens(e)\t(9)temp(e)\t(10)tprim(e)\t(11)fprim(e)\n')
##    for i in range(0, shape(profdata)[0]):
##        f.write(str(format2(s[i])))       
##        f.write(str(format2(nine[i])))
##        f.write(str(format2(tite[i])))
##        f.write(str(format2(dens_i[i])))
##        f.write(str(format2(t_i[i])))
##        f.write(str(format2(tprim_i[i])))
##        f.write(str(format2(fprim_i[i])))
##        f.write(str(format2(dens_e[i])))
##        f.write(str(format2(t_e[i])))
##        f.write(str(format2(tprim_e[i])))
##        f.write(str(format2(fprim_e[i])+'\n'))

