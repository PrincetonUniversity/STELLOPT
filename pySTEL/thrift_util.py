#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Main routine
if __name__=="__main__":
	import sys
	from argparse import ArgumentParser
	import numpy as np
	import matplotlib.pyplot as pyplot
	from libstell.thrift import THRIFT
	parser = ArgumentParser(description= 
		'''Provides a tool for making simple plots of THRIFT results.''')
	parser.add_argument("-t", "--thrift", dest="thrift_ext",
		help="THRIFT output HDF5 file.", default = None)
	parser.add_argument("-p", "--plot", dest="lplot", action='store_true',
		help="Plot the THRIFT file.", default = False)
	parser.add_argument("--print_bootstrap", dest="lprint_bootstrap", action='store_true',
		help="Print the boostrap current.", default = False)
	parser.add_argument("--print_potential", dest="lprint_potential", action='store_true',
		help="Print the electrostatic potential.", default = False)
	parser.add_argument("--print_profiles", dest="lprint_profiles", action='store_true',
		help="Print the plasma profiles.", default = False)
	parser.add_argument("--time", dest="timestamp",
		help="Provide a timestamp.", default = None, type=float)
	parser.add_argument("--save", dest="lsave", action='store_true',
		help="Save the plots with ext names.", default = False)
	args = parser.parse_args()
	thrift_data = THRIFT()
	loutput = False
	if args.thrift_ext:
		try:
			thrift_data.read_thrift(args.thrift_ext)
			loutput = True
			thrift_ext = args.thrift_ext.split('thrift_')[1]
			thrift_ext = thrift_ext.split('.h5')[0]
		except:
			print(f'Could not THRIFT output file: {args.thrift_ext}')
		# Make a general time plot
		if (loutput):
			if (args.lplot and (type(args.timestamp) is type(None))):
				color_order = ['b','r','g','c','m','y']
				px = 1/pyplot.rcParams['figure.dpi']
				fig=pyplot.figure(figsize=(1024*px,768*px))
				ax=fig.add_subplot(221)
				axb = ax.twinx()
				pyplot.subplots_adjust(hspace=0.3,wspace=0.4)
				nspecies = thrift_data.THRIFT_TEMP.shape[2]
				for i in range(nspecies):
					ax.plot(thrift_data.THRIFT_T,thrift_data.THRIFT_TEMP[:,0,i]/1000.,color=color_order[i])
					axb.plot(thrift_data.THRIFT_T,thrift_data.THRIFT_DENS[:,0,i]/1.0E20,'--',color=color_order[i])
				ax.legend()
				ax.set_xlabel('Time [s]')
				ax.set_ylabel('Temperature [keV]')
				ax.set_title('Profile Evolution')
				ax.set_ylim([0.0,25])
				ax.set_xlim([0,max(thrift_data.THRIFT_T)])
				axb.set_ylabel('Density $x10^{20}$ [$m^{-3}$]')
				axb.set_ylim([0.0,3])
				axb.set_xlim([0,max(thrift_data.THRIFT_T)])
				ax=fig.add_subplot(222)
				ax.plot(thrift_data.THRIFT_T,thrift_data.THRIFT_I[:,-1]/1000.,label='Total Current')
				ax.plot(thrift_data.THRIFT_T,thrift_data.THRIFT_IPLASMA[:,-1]/1000.,label='Plasma Current')
				ax.plot(thrift_data.THRIFT_T,thrift_data.THRIFT_ISOURCE[:,-1]/1000.,label='Source Current')
				ax.set_xlabel('Time [s]')
				ax.set_ylabel('Current [kA]')
				ax.set_title('Current Evolution')
				ax.set_xlim([0,max(thrift_data.THRIFT_T)])
				ax.legend()
				ax=fig.add_subplot(223)
				ax.plot(thrift_data.THRIFT_T,thrift_data.THRIFT_IOTA[:,0],label='Iota (core)')
				ax.plot(thrift_data.THRIFT_T,thrift_data.THRIFT_IOTA[:,-1],label='Iota (edge)')
				ax.set_xlabel('Time [s]')
				ax.set_ylabel('Rotational Transform')
				ax.set_title('Iota Evolution')
				ax.set_xlim([0,max(thrift_data.THRIFT_T)])
				ax.legend()
				ax=fig.add_subplot(224)
				ax.plot(thrift_data.THRIFT_T,thrift_data.THRIFT_IBOOT[:,-1]/1000.,label='Bootstrap Current')
				ax.plot(thrift_data.THRIFT_T,thrift_data.THRIFT_IECCD[:,-1]/1000.,label='ECCD')
				ax.plot(thrift_data.THRIFT_T,thrift_data.THRIFT_INBCD[:,-1]/1000.,label='NBI')
				ax.plot(thrift_data.THRIFT_T,thrift_data.THRIFT_IOHMIC[:,-1]/1000.,label='Ohmic')
				ax.set_xlabel('Time [s]')
				ax.set_ylabel('Current [kA]')
				ax.set_title('Current Source Evolution')
				ax.set_xlim([0,max(thrift_data.THRIFT_T)])
				ax.legend()
				pyplot.show()
				if (args.lsave): fig.savefig(f'thrift_evolution_{thrift_ext}.png', dpi=fig.dpi)
			elif (args.lplot):
				# For now we just plot the end conditino
				color_order = ['b','r','g','c','m','y']
				px = 1/pyplot.rcParams['figure.dpi']
				pyplot.rc('font', size=14)
				fig=pyplot.figure(figsize=(1024*px,768*px))
				ax=fig.add_subplot(221)
				axb = ax.twinx()
				pyplot.subplots_adjust(hspace=0.3,wspace=0.4)
				nspecies = thrift_data.THRIFT_TEMP.shape[2]
				THRIFT_RHO = np.sqrt(thrift_data.THRIFT_S)
				for i in range(nspecies):
					s,temp = thrift_data.get_temperature_prof(species=i,time=args.timestamp)
					rho = np.sqrt(s)
					ax.plot(rho,temp/1000.,color=color_order[i])
					s,dens = thrift_data.get_density_prof(species=i,time=args.timestamp)
					rho = np.sqrt(s)
					axb.plot(rho,dens/1.0E20,'--',color=color_order[i])
				ax.legend()
				ax.set_xlabel('Norm. Minor Raidus (r/a)')
				ax.set_ylabel('Temperature [keV]')
				ax.set_title('Kinetic Profiles')
				ax.set_ylim([0.0,25])
				ax.set_xlim([0,1.0])
				axb.set_ylabel('Density $x10^{20}$ [$m^{-3}$]')
				axb.set_ylim([0.0,3])
				axb.set_xlim([0,1.0])
				ax=fig.add_subplot(222)
				s_j,jtotal = thrift_data.get_j_prof(time=args.timestamp)
				s_jboot,jboot = thrift_data.get_jboot_prof(time=args.timestamp)
				s_jsource,jsource = thrift_data.get_jsource_prof(time=args.timestamp)
				ax.plot(np.sqrt(s_j),jtotal/1000.,label='Total Current')
				ax.plot(np.sqrt(s_jboot),jboot/1000.,label='Plasma Current')
				ax.plot(np.sqrt(s_jsource),jsource/1000.,label='Source Current')
				ax.set_xlabel('Norm. Minor Raidus (r/a)')
				ax.set_ylabel('Current Density [kA/m^2]')
				ax.set_title('Current Density Profile')
				ax.set_xlim([0,1.0])
				ax.text(0.02,0.12,rf'I_TOTAL: {thrift_data.get_I_total(time=args.timestamp):6.1f}', horizontalalignment='left',\
					verticalalignment='center', transform=ax.transAxes)
				ax.text(0.02,0.05,rf'I_BOOT: {thrift_data.get_Iboot_total(time=args.timestamp):6.1f}', horizontalalignment='left',\
					verticalalignment='center', transform=ax.transAxes)
				ax.legend()
				ax=fig.add_subplot(223)
				s,iota = thrift_data.get_iota_prof(time=args.timestamp)
				ax.plot(np.sqrt(s),iota,label='Iota (core)')
				ax.set_xlabel('Norm. Minor Raidus (r/a)')
				ax.set_ylabel('Rotational Transform')
				ax.set_title('Iota Profile')
				ax.set_xlim([0,1.0])
				ax.legend()
				ax=fig.add_subplot(224)
				s,er = thrift_data.get_er_prof(time=args.timestamp)
				ax.plot(np.sqrt(s),er/1000.,label=r'$E_R$')
				ax.set_xlabel('Norm. Minor Raidus (r/a)')
				ax.set_ylabel('Radial Electric Field [kV/m]')
				ax.set_title('Radial Electric Profile')
				ax.set_xlim([0,1.0])
				ax.legend()
				pyplot.show()
				if (args.lsave): fig.savefig(f'thrift_slice_{thrift_ext}_{int(args.timestamp*1000):0d}ms.png', dpi=fig.dpi)
			if (args.lprint_bootstrap and loutput):
				Ibs   = thrift_data.get_Iboot_total(time=args.timestamp) # in kA
				s,jbs = thrift_data.get_jboot_prof(time=args.timestamp)
				print(f"  CURTOR = {Ibs*1000.:20.10E}")
				temp = ' '.join(format(k, '8.5f') for k in s)
				print(f"  AC_AUX_S = "+temp)
				temp = ' '.join(format(k, '20.10E') for k in jbs)
				print(f"  AC_AUX_F = "+temp)
			if (args.lprint_potential and loutput):
				s,pot = thrift_data.get_pot_prof(time=args.timestamp)
				temp = ' '.join(format(k, '8.5f') for k in s)
				print(f"  POT_AUX_S = "+temp)
				temp = ' '.join(format(k, '20.10E') for k in pot)
				print(f"  POT_AUX_F = "+temp)
			if (args.lprint_profiles and loutput):
				s,ne = thrift_data.get_density_prof(species=0,time=args.timestamp)
				temp = ' '.join(format(k, '8.5f') for k in s)
				print(f"  NE_AUX_S = "+temp)
				temp = ' '.join(format(k, '20.10E') for k in ne)
				print(f"  NE_AUX_F = "+temp)
				s,te = thrift_data.get_temperature_prof(species=0,time=args.timestamp)
				temp = ' '.join(format(k, '8.5f') for k in s)
				print(f"  TE_AUX_S = "+temp)
				temp = ' '.join(format(k, '20.10E') for k in te)
				print(f"  TE_AUX_F = "+temp)
	sys.exit(0)