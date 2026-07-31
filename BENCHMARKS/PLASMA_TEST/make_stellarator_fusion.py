#!/usr/bin/env python3
import numpy as np
import sys
import joblib
from types import SimpleNamespace
from scipy.interpolate import interp1d
from libstell.vmec import VMEC
from libstell.plasma_solver import PLASMA_SOLVER
from libstell.plasma import PLASMA

plasma = PLASMA(['electrons','deuterium','tritium','helium4'])

solver = PLASMA_SOLVER(['electrons','deuterium','tritium','helium4'],tau_fast_alphas=0.1,constrain_nT=True,solve_fast_alphas=True)

#################################################################################################
######################### set initial profiles ##################################################
rho_init = np.linspace(0,1,100)

# constant densities in the beggining
NE0_init = 5E19 # can't be smaller than 10eV otherise clog gets negative
NE_edge_init = NE0_init
NE_init = NE_edge_init + (NE0_init-NE_edge_init)*(1-rho_init**2)

solver.set_initial_profile('density', 'electrons', rho_vals=rho_init, profile_vals=NE_init)
solver.set_initial_profile('density', 'deuterium', rho_vals=rho_init, profile_vals=0.495*NE_init)
solver.set_initial_profile('density', 'tritium',   rho_vals=rho_init, profile_vals=0.495*NE_init)
solver.set_initial_profile('density', 'helium4',   rho_vals=rho_init, profile_vals=0.005*NE_init)

# constant temperatures in the beggining
TE0_init = 200
TE_edge_init = TE0_init
TE_init = TE_edge_init + (TE0_init-TE_edge_init)*(1-rho_init**2)

solver.set_initial_profile('temperature', 'electrons', rho_vals=rho_init, profile_vals=TE_init)
solver.set_initial_profile('temperature', 'deuterium', rho_vals=rho_init, profile_vals=TE_init)
solver.set_initial_profile('temperature', 'tritium',   rho_vals=rho_init, profile_vals=TE_init)
solver.set_initial_profile('temperature', 'helium4',   rho_vals=rho_init, profile_vals=TE_init)

#################################################################################################
###################### set Boundary Conditions ##################################################

solver.set_edge_boundary_condition('density', 'electrons',NE_init[-1])
solver.set_edge_boundary_condition('density', 'deuterium',0.495*NE_init[-1])
solver.set_edge_boundary_condition('density', 'tritium'  ,0.495*NE_init[-1])
solver.set_edge_boundary_condition('density', 'helium4'  ,0.005*NE_init[-1])

solver.set_edge_boundary_condition('temperature', 'electrons', TE_init[-1])
solver.set_edge_boundary_condition('temperature', 'deuterium', TE_init[-1])
solver.set_edge_boundary_condition('temperature', 'tritium'  , TE_init[-1])
solver.set_edge_boundary_condition('temperature', 'helium4'  , TE_init[-1])

#################################################################################################
###################### set Heating and Fueling Scenarios#########################################
time_vec  = [0.00,5.00,10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0,50.0,55.0,60.0,65.0,70.0,75.0,80.0,85.0,90.0,95.0,100.0,130.0,150.0,160.0,170.0,180.0,190.0,200.0]
ecrh_vec  = [5.00,5.00,10.0,10.0,10.0,10.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0,50.0,60.0,70.0,80.0,80.0,80.0,80.00,80.00,80.00,60.00,40.00,20.00,0.000,0.000]
fuel_vec  = [5.00,5.00,5.00,5.00,5.00,10.0,10.0,10.0,10.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,22.0,22.00,22.00,22.00,22.00,22.00,22.00,22.00,22.00]

fact_ECRH = interp1d(time_vec,ecrh_vec,
	                  kind='previous', bounds_error=False, fill_value=(0.0,0.0))

fact_fuel = interp1d(time_vec,np.array(fuel_vec)*1.0E19,
	            kind='linear', bounds_error=False, fill_value=(0.0,0.0))


#################################################################################################
###################### set Energy Sources #######################################################
solver.set_energy_source('electrons','time_dependent_gaussian', total_power=1E6, sigma_rho=0.10, rho_0=0.0, time_dependent_factor=fact_ECRH)

solver.set_energy_source('electrons','Coll_Heat_Exchange')
solver.set_energy_source('deuterium','Coll_Heat_Exchange')
solver.set_energy_source('tritium','Coll_Heat_Exchange')
solver.set_energy_source('helium4','Coll_Heat_Exchange')

solver.set_energy_source('electrons','alpha_heating',fraction_alpha_heating=0.85)
solver.set_energy_source('deuterium','alpha_heating',fraction_alpha_heating=0.075)
solver.set_energy_source('tritium','alpha_heating',fraction_alpha_heating=0.075)

solver.set_energy_source('electrons','Bremsstrahlung')

#################################################################################################
###################### set Particle Sources #######################################################


#################################################################################################
###################### set Particle Sources #######################################################
K_fact = 50.0*20.0
I_fact = (1.2/2)*10.0
D_fact = 1.0

# This is an edge source to mimimc gas fueling (these must be the same)
solver.set_particle_source('deuterium','PID_edense_gaussian', 
	rho_0=1.00, max_injected_particles_per_sec = 1.0E24, sigma_rho = 0.025,
	pidK=K_fact,pidI=I_fact,pidD=D_fact,noise_level=0.00,
	time_dependent_electron_dens_axis=fact_fuel)
solver.set_particle_source('tritium','PID_edense_gaussian', 
	rho_0=1.00, max_injected_particles_per_sec = 1.0E24, sigma_rho = 0.025,
	pidK=K_fact,pidI=I_fact,pidD=D_fact,noise_level=0.00,
	time_dependent_electron_dens_axis=fact_fuel)

solver.set_particle_source('helium4'  ,'fast_alphas_source')
solver.set_particle_source('deuterium','alpha_particles_sink')
solver.set_particle_source('tritium','alpha_particles_sink')

#################################################################################################
###################### set Equilibrium ##########################################################

solver.set_equilibrium(type='cylindrical', aminor=2.0, Rmajor=20.0, B=6.0, iota23=0.9)

#################################################################################################
###################### set Heat Fluxes ##########################################################

solver.set_heat_fluxes(type='beurskens', chi_base=0.03, aLT_critical=2.20, alpha=1.0, stiffness=0.4, chi_electrons=0.5, convective_fact=0.0,mass_ref_species=(plasma.mass_database['deuterium']+plasma.mass_database['tritium'])/2.0)

####################################################################################################
###################### set Particle Fluxes ##########################################################
Dn = 0.1
solver.set_particle_fluxes(type='diffusive',Dn=Dn)

####################################################################################################
####################################  RUN ##########################################################

DT=5E-3
solver.run(Nr=101,dt=DT,tstart=0.0,tend=200.0,tolerance=1E-3,max_subiter=20,output_filename='stellarator_fusion')
