
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec 
from stellapy.data.write_profile import calculate_gradients, guess_columns
from stellapy.utils import get_filesInFolder
from stellapy.plot.utils import load_plotbox2d, load_styleFigures

def plot_profilesVsRho(\
            # Get the data 
            raw_files, \
            # Define the radial coordinate
            x1='rho', \
            a1=None, \
            x2='rho', \
            a2=None, \
            # State in which columns the data can be found
            x_col1=None, \
            dens_col1=None, \
            Ti_col1=None, \
            Te_col1=None, \
            x_col2=None, \
            dens_col2=None, \
            Ti_col2=None, \
            Te_col2=None, \
            # If it is called from the GUI, the figure exists
            save = False, \
            ax1 = None, ax1_twin = None, \
            ax2 = None, ax2_twin = None, \
            ax3 = None, ax3_twin = None, \
            plot = "profiles"):
    ''' Obtain the gradients and values of the n_e, T_e and T_i profiles at specific rho values 

    Parameters
    ----------
    folder, source_file : str
        specify the folder in the runs directory and the source_file with extension
    x : {'rho', 's', 'r'}
        Specify the units of the x-axis, this will be converted into rho
        If [x] = 'r' also specify the minor radius [a]
    x_col, dens_col, Ti_col, Te_col : int
        Specify in which columns the quantities can be found
    plot : {"a/Ln", "profiles", "gradients", "normalized gradients", "relative comparison"}
        If the figure already exists, select what to plot on it

    BASH: plot_profiles --x=0 --n=2 --Te=3 --Ti=4
    '''    
        
    # Plot everything if the figure doesn't exist
    if ax1 == None:
        plot = "all"
        plt.rc('font', size=25)
    
    # Set-up the figure for the profiles
    if plot in ["profiles", "gradients", "normalized gradients"]:
        # Profiles of density; electron temperature; ion temperature
        if plot in ["profiles"]: 
            ax_profileDensity = ax1 
            ax_profileEleTemp = ax2
            ax_profileIonTemp = ax3
        # Real gradients of density; electron temperature; ion temperature
        if plot in ["gradients"]: 
            ax_gradientsDensity = ax1 
            ax_gradientsEleTemp = ax2
            ax_gradientsIonTemp = ax3 
        # Normalized logaritmic gradients of density; electron temperature; ion temperature
        if plot in ["normalized gradients"]: 
            ax_normGradDensity = ax1
            ax_normGradEleTemp = ax2
            ax_normGradIonTemp = ax3
    if plot == "all":
        fig = plt.figure(figsize=(13, 10))
        gs = gridspec.GridSpec(3, 3)
        gs.update(wspace=0.35, hspace=0.15, top=0.97, bottom=0.1)
        ax_profileDensity = plt.subplot(gs[0])
        ax_profileEleTemp = plt.subplot(gs[1])
        ax_profileIonTemp = plt.subplot(gs[2])
        ax_gradientsDensity = plt.subplot(gs[3])
        ax_gradientsEleTemp = plt.subplot(gs[4])
        ax_gradientsIonTemp = plt.subplot(gs[5])
        ax_normGradDensity = plt.subplot(gs[6])
        ax_normGradEleTemp = plt.subplot(gs[7])
        ax_normGradIonTemp = plt.subplot(gs[8])
        
    # Make the plots pretty
    if plot=="all" or plot in ["profiles"]:
        xlab = "" if plot=="all" else  '$\\rho = r/a$'
        load_plotbox2d(x_range=[0,1], x_label=xlab, y_label='${n_e}$ $[10^{19}$m$^{-3}]$', title="", ax=ax_profileDensity, title_size='small')
        load_plotbox2d(x_range=[0,1], x_label=xlab, y_label='${T_i}$ [keV]', title="", ax=ax_profileEleTemp, title_size='small')
        load_plotbox2d(x_range=[0,1], x_label=xlab, y_label='${T_e}$ [keV]', title="", ax=ax_profileIonTemp, title_size='small')
    if plot=="all" or plot in ["gradients"]:
        xlab = "" if plot=="all" else  '$\\rho = r/a$'    
        load_plotbox2d(x_range=[0,1], x_label=xlab, y_label="$-n'_e$ $[10^{19}$m$^{-4}]$", ax=ax_gradientsDensity)
        load_plotbox2d(x_range=[0,1], x_label=xlab, y_label="$-T'_i$ [keV/m]", ax=ax_gradientsEleTemp)
        load_plotbox2d(x_range=[0,1], x_label=xlab, y_label="$-T'_e$ [keV/m]", ax=ax_gradientsIonTemp)
    if plot=="all" or plot in ["normalized gradients"]:
        xlab = '$\\rho = r/a$'
        load_plotbox2d(x_range=[0,1], x_label=xlab, y_label='$a/L_{n_e}$', ax=ax_normGradDensity)
        load_plotbox2d(x_range=[0,1], x_label=xlab, y_label='$a/L_{T_i}$', ax=ax_normGradEleTemp)
        load_plotbox2d(x_range=[0,1], x_label=xlab, y_label='$a/L_{T_e}$', ax=ax_normGradIonTemp)
        
    # All normalized gradients
    if plot in ["a/Ln"]:
        load_styleFigures()
        ax_allNormGradProf1 = ax1
        ax_allNormGradProf2 = ax2
        load_plotbox2d(x_label='$\\rho = r/a$', y_label='$a/L$', fig_size=(8.5, 7.5), ax=ax_allNormGradProf1)
        load_plotbox2d(x_label='$\\rho = r/a$', y_label='$a/L$', fig_size=(8.5, 7.5), ax=ax_allNormGradProf2)
        ax_allNormGradProf1Twin = ax1_twin    # Create a second axis for Te/Ti
        ax_allNormGradProf2Twin = ax2_twin    # Create a second axis for Te/Ti
        ax_allNormGradProf1Twin.set_ylabel("$T_e/T_i$", color='navy'); ax_allNormGradProf1Twin.tick_params(axis='y', labelcolor='navy', colors='navy')
        ax_allNormGradProf2Twin.set_ylabel("$T_e/T_i$", color='navy'); ax_allNormGradProf2Twin.tick_params(axis='y', labelcolor='navy', colors='navy')
        ax_allNormGradProf = ax3
        load_plotbox2d(x_label='$\\rho = r/a$', y_label='$a/L$', fig_size=(12, 6), ax=ax_allNormGradProf)
        ax_allNormGradProfTwin = ax3_twin 
        ax_allNormGradProfTwin.set_ylabel("$T_e/T_i$", color='navy'); ax_allNormGradProfTwin.tick_params(axis='y', labelcolor='navy')
        
    # Compare two profiles
    if plot in ["relative comparison"]:
        ax_allNormGradProf = ax1
        load_plotbox2d(x_label='$\\rho = r/a$', y_label='$a/L$', fig_size=(12, 6), ax=ax_allNormGradProf)
        ax_allNormGradProfTwin = ax1_twin   
        ax_allNormGradProfTwin.set_ylabel("$T_e/T_i$", color='navy'); ax_allNormGradProfTwin.tick_params(axis='y', labelcolor='navy')
        
        ax_relativeNormGrad = ax2
        load_plotbox2d(x_label='$\\rho = r/a$', y_label='$( a/L_{2} - a/L_{1} ) / (a/L_{1})$', ax=ax_relativeNormGrad)
        ax_relativeNormGradTwin = ax2_twin   
        ax_relativeNormGradTwin.set_ylabel("$( (T_e/T_i)_{2} - (T_e/T_i)_{1} ) / (T_e/T_i)_{1}$", color='navy'); ax_relativeNormGradTwin.tick_params(axis='y', labelcolor='navy')

    if plot in ["all"]:
        fig3 = plt.figure(figsize=(7, 6))
        gs2  = gridspec.GridSpec(1, 2)
        ax_allNormGradProf1 = plt.subplot(gs2[0]) 
        ax_allNormGradProf2 = plt.subplot(gs2[1]) 
        load_plotbox2d(x_label='$\\rho = r/a$', y_label='$a/L$', fig_size=(8.5, 7.5), ax=ax_allNormGradProf1)
        load_plotbox2d(x_label='$\\rho = r/a$', y_label='$a/L$', fig_size=(8.5, 7.5), ax=ax_allNormGradProf2)
        ax_allNormGradProf1Twin = ax_allNormGradProf1.twinx()    # Create a second axis for Te/Ti
        ax_allNormGradProf2Twin = ax_allNormGradProf2.twinx()    # Create a second axis for Te/Ti
        ax_allNormGradProf1Twin.set_ylabel("$T_e/T_i$", color='navy'); ax_allNormGradProf1Twin.tick_params(axis='y', labelcolor='navy')
        ax_allNormGradProf2Twin.set_ylabel("$T_e/T_i$", color='navy'); ax_allNormGradProf2Twin.tick_params(axis='y', labelcolor='navy')
        
        fig2 = plt.figure(figsize=(14, 6))
        gs3 = gridspec.GridSpec(1, 1)
        ax_allNormGradProf = plt.subplot(gs3[0])  
        load_plotbox2d(x_label='$\\rho = r/a$', y_label='$a/L$', fig_size=(12, 6), ax=ax_allNormGradProf)
        ax_allNormGradProfTwin = ax_allNormGradProf.twinx()    
        ax_allNormGradProfTwin.set_ylabel("$T_e/T_i$", color='navy'); ax_allNormGradProfTwin.tick_params(axis='y', labelcolor='navy')
        
        fig4 = plt.figure(figsize=(7.5, 6))
        gs4  = gridspec.GridSpec(1, 1)
        ax_relativeNormGrad = plt.subplot(gs4[0])  
        load_plotbox2d(x_label='$\\rho = r/a$', y_label='$( a/L_{2} - a/L_{1} ) / (a/L_{1})$', ax=ax_relativeNormGrad)
        ax_relativeNormGradTwin = ax_relativeNormGrad.twinx()    
        ax_relativeNormGradTwin.set_ylabel("$( (T_e/T_i)_{2} - (T_e/T_i)_{1} ) / (T_e/T_i)_{1}$", color='navy'); ax_relativeNormGradTwin.tick_params(axis='y', labelcolor='navy')

    # Count the data files
    count = 0

    # Read the data file containig the profiles
    for raw_file in raw_files:
        
        # Count the data files
        count += 1
        
        # Read the data
        try: data = np.loadtxt(raw_file, dtype='float')
        except:
            data = np.loadtxt(raw_file, dtype=np.str, delimiter='\t', skiprows=1)
            data = np.char.replace(data, ',', '.').astype(np.float64)

        # If the columns are not specified, try to guess
        if count==1:
            x = x1; a = a1
            if (x_col1 is None) or (dens_col1 is None) or (Ti_col1 is None) or (Te_col1 is None):
                x_col, dens_col, Ti_col, Te_col  = guess_columns(raw_file)
            else:
                x_col, dens_col, Ti_col, Te_col = x_col1, dens_col1, Ti_col1, Te_col1
        if count==2:
            x = x2; a = a2
            if (x_col2 is None) or (dens_col2 is None) or (Ti_col2 is None) or (Te_col2 is None):
                x_col, dens_col, Ti_col, Te_col = guess_columns(raw_file)
            else:
                x_col, dens_col, Ti_col, Te_col = x_col2, dens_col2, Ti_col2, Te_col2
            
        # The "Ne_fits" file contains the data of two experiments!
        if "Ne_fits" in str(raw_file):
            if count==1: x_col=0; dens_col=1; Ti_col=3; Te_col=2
            if count==2: x_col=0; dens_col=4; Ti_col=6; Te_col=5
        
        # Change the unit along the x-axis to rho=r/a=sqrt(s)
        if x == 'r':  data[:,x_col] = data[:,x_col]/a
        if x == 's':  data[:,x_col] = np.sqrt(data[:,x_col])

        # Add shot labels
        if count==1:
            label_S = "Profile 1"
            color   = "navy"
            style   = '-'
        if count==2:
            label_S = "Profile 2"  
            color   = "crimson"
            style   = '--'
                
        # Density: Get the profiles and gradients and plot them
        Y_prof, Y_grad, a_over_Ly, x_vec = calculate_gradients(data,rho_col=x_col,Y_col=dens_col)
        if plot=="all" or plot in ["profiles"]:
            ax_profileDensity.plot(x_vec, Y_prof, style, color=color, linewidth=4, markerfacecolor='white', label=label_S)
        if plot=="all" or plot in ["gradients"]:
            ax_gradientsDensity.plot(x_vec, -Y_grad, style, color=color, linewidth=4, markerfacecolor='white', label=label_S)
        if plot=="all" or plot in ["normalized gradients"]:    
            ax_normGradDensity.plot(x_vec, a_over_Ly, style, color=color, linewidth=4, markerfacecolor='white', label=label_S)
        if count==1: a_over_Ln_ld = a_over_Ly; n_prof_ld = Y_prof
        if count==2: a_over_Ln_hd = a_over_Ly; n_prof_hd = Y_prof
 
        # Ti: Get the profiles and gradients and plot them
        Y_prof, Y_grad, a_over_Ly, x_vec = calculate_gradients(data,rho_col=x_col,Y_col=Ti_col)
        if plot=="all" or plot in ["profiles"]:
            ax_profileEleTemp.plot(x_vec, Y_prof,    style, color=color, linewidth=4, markerfacecolor='white', label=label_S)
        if plot=="all" or plot in ["gradients"]:
            ax_gradientsEleTemp.plot(x_vec, -Y_grad,   style, color=color, linewidth=4, markerfacecolor='white', label=label_S)
        if plot=="all" or plot in ["normalized gradients"]:
            ax_normGradEleTemp.plot(x_vec, a_over_Ly, style, color=color, linewidth=4, markerfacecolor='white', label=label_S)
        if count==1: a_over_LTi_ld = a_over_Ly; Ti_prof_ld = Y_prof
        if count==2: a_over_LTi_hd = a_over_Ly; Ti_prof_hd = Y_prof        

        # Te: Get the profiles and gradients and plot them
        Y_prof, Y_grad, a_over_Ly, x_vec = calculate_gradients(data,rho_col=x_col,Y_col=Te_col)
        if plot=="all" or plot in ["profiles"]:
            ax_profileIonTemp.plot(x_vec, Y_prof,    style, color=color, linewidth=4, markerfacecolor='white', label=label_S)
        if plot=="all" or plot in ["gradients"]:
            ax_gradientsIonTemp.plot(x_vec, -Y_grad,   style, color=color, linewidth=4, markerfacecolor='white', label=label_S)
        if plot=="all" or plot in ["normalized gradients"]:
            ax_normGradIonTemp.plot(x_vec, a_over_Ly, style, color=color, linewidth=4, markerfacecolor='white', label=label_S)
        if count==1: a_over_LTe_ld = a_over_Ly; Te_prof_ld = Y_prof
        if count==2: a_over_LTe_hd = a_over_Ly; Te_prof_hd = Y_prof     
       
    # Plot all gradients together in two subplots, one shot per subplot
    if plot=="all" or plot in ["a/Ln"]:
        ax_allNormGradProf1.plot(x_vec, a_over_Ln_ld, color='black', linewidth=3, markerfacecolor='white', label="$a/L_n$")
        ax_allNormGradProf1.plot(x_vec, a_over_LTi_ld, color='green', linewidth=3, markerfacecolor='white', label="$a/L_{T_i}$")
        ax_allNormGradProf1.plot(x_vec, a_over_LTe_ld, color='red', linewidth=3, markerfacecolor='white', label="$a/L_{T_e}$")
        ax_allNormGradProf1Twin.plot(x_vec, Te_prof_ld/Ti_prof_ld, color='navy', linewidth=3, markerfacecolor='white', label='')
        ax_allNormGradProf1.plot(-1, -1, color='navy', linewidth=3, label='$T_e/T_i$')
        if len(raw_files)!=1: 
            ax_allNormGradProf2.plot(x_vec, a_over_Ln_hd, color='black', linewidth=3, markerfacecolor='white', label="$a/L_n$")
            ax_allNormGradProf2.plot(x_vec, a_over_LTi_hd, color='green', linewidth=3, markerfacecolor='white', label="$a/L_{T_i}$")
            ax_allNormGradProf2.plot(x_vec, a_over_LTe_hd, color='red', linewidth=3, markerfacecolor='white', label="$a/L_{T_e}$")
            ax_allNormGradProf2Twin.plot(x_vec, Te_prof_hd/Ti_prof_hd, color='navy', linewidth=3, markerfacecolor='white', label='')
            ax_allNormGradProf2.plot(-1, -1, color='navy', linewidth=3, label='$T_e/T_i$')

    # Plot al gradients together in one subplot
    if plot=="all" or plot in ["relative comparison", "a/Ln"]:
        ax_allNormGradProf.plot(x_vec, a_over_Ln_ld, color='black', ls='-', linewidth=3, markerfacecolor='white', label="$a/L_n$")
        ax_allNormGradProf.plot(x_vec, a_over_LTi_ld, color='green', ls='-', linewidth=3, markerfacecolor='white', label="$a/L_{T_i}$")
        ax_allNormGradProf.plot(x_vec, a_over_LTe_ld, color='red', ls='-', linewidth=3, markerfacecolor='white', label="$a/L_{T_e}$")
        ax_allNormGradProfTwin.plot(x_vec, Te_prof_ld/Ti_prof_ld, color='navy', ls='-', linewidth=3, markerfacecolor='white', label='')
        ax_allNormGradProf.plot(-1, -1, color='navy', linewidth=3, label='$T_e/T_i$')
        if len(raw_files)!=1: 
            ax_allNormGradProf.plot(x_vec, a_over_Ln_hd, color='black', ls='--', linewidth=3, markerfacecolor='white', label="")
            ax_allNormGradProf.plot(x_vec, a_over_LTi_hd, color='green', ls='--', linewidth=3, markerfacecolor='white', label="")
            ax_allNormGradProf.plot(x_vec, a_over_LTe_hd, color='red', ls='--', linewidth=3, markerfacecolor='white', label="")
            ax_allNormGradProfTwin.plot(x_vec, Te_prof_hd/Ti_prof_hd, color='navy', ls='--', linewidth=3, markerfacecolor='white', label='')
            ax_allNormGradProf.plot(-1, -1, color='navy', linewidth=3, label='')
    
    # Plot the percentual differences: Profile 2 is reference so difference = (ld-hd)/hd
    if plot=="all" or plot in ["relative comparison"]:
        if len(raw_files)!=1:
            ax_relativeNormGrad.plot(x_vec, (a_over_Ln_ld-a_over_Ln_hd)/a_over_Ln_hd*100, color='black', ls='-', linewidth=3, label="$a/L_n$") 
            ax_relativeNormGrad.plot(x_vec, (a_over_LTi_ld-a_over_LTi_hd)/a_over_LTi_hd*100, color='green', ls='-', linewidth=3, label="$a/L_{T_i}$") 
            ax_relativeNormGrad.plot(x_vec, (a_over_LTe_ld-a_over_LTe_hd)/a_over_LTe_hd*100, color='red', ls='-', linewidth=3, label="$a/L_{T_e}$") 
            ax_relativeNormGradTwin.plot(x_vec, (Te_prof_ld/Ti_prof_ld-Te_prof_hd/Ti_prof_hd)/(Te_prof_hd/Ti_prof_hd)*100, color='navy', ls='-', linewidth=3, label='')
            ax_relativeNormGrad.plot(-1, -1, color='navy', linewidth=3, label='$T_e/T_i$') 

    # Finish the figures
    if plot=="all" or plot in ["profiles"]:
        ax_profileDensity.set_title("Density profile")
        ax_profileEleTemp.set_title("Electron temperature profile")
        ax_profileIonTemp.set_title("Ion temperature profile")
        for ax in [ax_profileDensity, ax_profileEleTemp, ax_profileIonTemp]:
            ax.autoscale()
            ax.set_xlim(xmin=0, xmax=1)    
            ax.set_ylim(ymin=0)
            ax.set_ylim(ymax=10)
    if plot=="all" or plot in ["gradients"]:
        ax_gradientsDensity.set_title("Density gradients")
        ax_gradientsEleTemp.set_title("Electron temperature gradients")
        ax_gradientsIonTemp.set_title("Ion temperature gradients")
        for ax in [ax_gradientsDensity, ax_gradientsEleTemp, ax_gradientsIonTemp]:
            ax.autoscale()
            ax.set_xlim(xmin=0, xmax=1)    
            ax.set_ylim(ymin=0)
            ax.set_ylim(ymax=10)
    if plot=="all" or plot in ["normalized gradients"]:
        ax_normGradDensity.set_title("Normalized density gradients")
        ax_normGradEleTemp.set_title("Normalized electron temperature gradients")
        ax_normGradIonTemp.set_title("Normalized ion temperature gradients")
        for ax in [ax_normGradDensity, ax_normGradEleTemp, ax_normGradIonTemp]:
            ax.autoscale()
            ax.set_xlim(xmin=0, xmax=1)    
            ax.set_ylim(ymin=0)
            ax.set_ylim(ymax=10)
    if plot=="all" or plot in ["a/Ln"]:
        ax_allNormGradProf1.set_title("Profile 1")
        ax_allNormGradProf2.set_title("Profile 2")
        for ax in [ax_allNormGradProf1, ax_allNormGradProf2, ax_allNormGradProf1Twin, ax_allNormGradProf2Twin]:
            ax.autoscale()
            ax.set_xlim(xmin=0, xmax=1)        
            ax.set_ylim(ymin=0)
            ax.set_ylim(ymax=10)
    if plot=="all" or plot in ["relative comparison", "a/Ln"]:
        ax_allNormGradProf.set_title("Both profiles")
        for ax in [ax_allNormGradProf, ax_allNormGradProfTwin]:
            ax.autoscale()
            ax.set_xlim(xmin=0, xmax=1)        
            ax.set_ylim(ymin=0)
            ax.set_ylim(ymax=10)
    if plot=="all" or plot in ["relative comparison"]:   
        ax_relativeNormGrad.set_title("Relative comparison") 
        ax_relativeNormGrad.yaxis.set_major_formatter(ticker.PercentFormatter())
        ax_relativeNormGradTwin.yaxis.set_major_formatter(ticker.PercentFormatter())
        for ax in [ax_relativeNormGrad, ax_relativeNormGradTwin]:
            ax.autoscale()
            ax.set_xlim(xmin=0, xmax=1)   
            ax.set_ylim(ymin=-20, ymax=20)   
    
    # Legend
    if plot=="all" or plot in ["profiles", "gradients", "normalized gradients"]:
        if plot=="all" or plot in ["profiles"]:
            for ax in [ax_profileDensity, ax_profileEleTemp, ax_profileIonTemp]:
                if len(raw_files)>1 or plot=="all":
                    ax.legend(loc='upper right', labelspacing=0.1, prop={'size':20})
        if plot=="all" or plot in ["gradients"]:
            for ax in [ax_gradientsDensity, ax_gradientsEleTemp, ax_gradientsIonTemp]:
                if len(raw_files)>1 or plot=="all":
                    ax.legend(loc='upper left', labelspacing=0.1, prop={'size':20})
        if plot=="all" or plot in ["normalized gradients"]:
            for ax in [ax_normGradDensity, ax_normGradEleTemp, ax_normGradIonTemp]:
                if len(raw_files)>1 or plot=="all":
                    ax.legend(loc='upper left', labelspacing=0.1, prop={'size':20})
        if plot=="all":
            ax.yaxis.set_label_coords(-0.13,0.5) 
                
    if plot=="all" or plot in ["a/Ln"]:
        for ax in [ax_allNormGradProf1, ax_allNormGradProf2, ax_allNormGradProf]:
            ax.legend(loc='upper left', labelspacing=0.1, prop={'size':20})
        for ax in [ax_allNormGradProfTwin]:
            ax.plot([-1,-2],[-1,-2],color="black",ls='-',label="Profile 1")
            ax.plot([-1,-2],[-1,-2],color="black",ls='--',label="Profile 2")
            ax.legend(loc='upper right', labelspacing=0.1, prop={'size':20})
            
    if plot=="all" or plot in ["relative comparison"]:
        for ax in [ax_relativeNormGrad]:
            ax.legend(labelspacing=0.1, prop={'size':20})           

    # Show the figure
    if plot=="all": plt.show()
    return 



