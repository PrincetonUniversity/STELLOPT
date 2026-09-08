#!/usr/bin/env python3
import numpy as np
import sys
import joblib
from types import SimpleNamespace
from scipy.interpolate import interp1d
from libstell.vmec import VMEC
from libstell.plasma_solver import PLASMA_SOLVER
from libstell.plasma import PLASMA

plasma = PLASMA(['electrons','hydrogen'])

solver = PLASMA_SOLVER(['electrons','hydrogen'],tau_fast_alphas=0.0,constrain_nT=False,solve_fast_alphas=False)

#################################################################################################
######################### set initial profiles ##################################################
rho_init = np.linspace(0,1,100)

# constant densities in the beggining
NE0_init = 1.1E19 # can't be smaller than 10eV otherise clog gets negative
NE_edge_init = 5E18
NE_init = NE_edge_init + (NE0_init-NE_edge_init)*(1-rho_init**6)

solver.set_initial_profile('density', 'electrons', rho_vals=rho_init, profile_vals=NE_init)
solver.set_initial_profile('density', 'hydrogen', rho_vals=rho_init, profile_vals=NE_init)

# constant temperatures in the beggining
TE0_init = 3000 # can't be smaller than 10eV otherise clog gets negative
TE_edge_init = 50
TE_init = TE_edge_init + (TE0_init-TE_edge_init)*(1-rho_init**2)

solver.set_initial_profile('temperature', 'electrons', rho_vals=rho_init, profile_vals=TE_init)
solver.set_initial_profile('temperature', 'hydrogen', rho_vals=rho_init, profile_vals=TE_init)

#################################################################################################
###################### set Boundary Conditions ##################################################

solver.set_edge_boundary_condition('density', 'electrons',NE_init[-1])
solver.set_edge_boundary_condition('density', 'hydrogen',NE_init[-1])

solver.set_edge_boundary_condition('temperature', 'electrons', TE_init[-1])
solver.set_edge_boundary_condition('temperature', 'hydrogen', TE_init[-1])

#################################################################################################
###################### set Heating and Fueling Scenarios#########################################
time_vec    = [0.00,1.00,1.10,1.60,1.70,3.50]
ecrh_vec    = [2.60,2.60,2.60,4.90,4.90,0.00]
fuel_vec    = [1.00,1.00,1.00,1.00,5.00,0.00]
pellet_vec  = [0.00,0.00,1.00,1.00,0.00,0.00]

fact_ECRH = interp1d(time_vec,ecrh_vec,
	                  kind='previous', bounds_error=False, fill_value=(0.0,0.0))

fact_fuel = interp1d(time_vec,np.array(fuel_vec)*1.0E19,
	            kind='linear', bounds_error=False, fill_value=(0.0,0.0))

fact_pellet = interp1d(time_vec,pellet_vec,
	                  kind='previous', bounds_error=False, fill_value=(0.0,0.0))


#################################################################################################
###################### set Energy Sources #######################################################
solver.set_energy_source('electrons','time_dependent_gaussian', total_power=1E6, sigma_rho=0.20, rho_0=0.0, time_dependent_factor=fact_ECRH)

solver.set_energy_source('electrons','Coll_Heat_Exchange')
solver.set_energy_source('hydrogen','Coll_Heat_Exchange')
solver.set_energy_source('electrons','Bremsstrahlung')

#################################################################################################
###################### set Particle Sources #######################################################
K_fact = 50.0*20.0
I_fact = (1.2/2)*10.0
D_fact = 1.0
# This is an edge source to mimimc gas fueling
solver.set_particle_source('hydrogen','PID_edense_gaussian', 
	rho_0=1.00, max_injected_particles_per_sec = 1.0E23, sigma_rho = 0.025,
	pidK=K_fact,pidI=I_fact,pidD=D_fact,noise_level=0.00,
	time_dependent_electron_dens_axis=fact_fuel)

# This is the pellet fueling
solver.set_particle_source('hydrogen','pellet_model', 
	pellet_size_mm=1.0, pellet_vel_ms=250.0, pellet_freq_Hz=15.0, 
	pellet_LtoD=1.0, pellet_mass_amu=1.0, pellet_density_kgm3=87.0,
	time_dependent_factor=fact_pellet)

#################################################################################################
###################### set Equilibrium ##########################################################

wout_file = 'wout_w7x.nc'
solver.set_equilibrium(type='VMEC', wout_path=wout_file)

#################################################################################################
###################### set Heat Fluxes ##########################################################

solver.set_heat_fluxes(type='beurskens', chi_base=0.03, aLT_critical=1.19, alpha=1.0, stiffness=0.89, chi_electrons=3.0, convective_fact=0.0,mass_ref_species=plasma.mass_database['hydrogen'])

####################################################################################################
###################### set Particle Fluxes ##########################################################
Dn = 0.09
solver.set_particle_fluxes(type='diffusive',Dn=Dn)

####################################################################################################
####################################  RUN ##########################################################

DT=1E-3
solver.run(Nr=101,dt=DT,tstart=0.0,tend=4.0,tolerance=1E-3,max_subiter=20,output_filename='w7x_h2_pellet')
