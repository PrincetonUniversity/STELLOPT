#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Main routine
if __name__=="__main__":
	import sys
	from argparse import ArgumentParser
	import numpy as np
	import matplotlib.pyplot as pyplot
	from libstell.vmec import VMEC, VMEC_INDATA
	from libstell.libstell import LIBSTELL
	parser = ArgumentParser(description= 
		'''Provides class for accessing vmec data also serves as a
		   simple tool for assessing vmec wout or input files.''')
	parser.add_argument("-v", "--vmec", dest="vmec_ext",
		help="VMEC file extension", default = None)
	parser.add_argument("-p", "--plot", dest="lplot", action='store_true',
		help="Plot the VMEC file.", default = False)
	parser.add_argument("-b", "--boozer", dest="lbooz", action='store_true',
		help="Output the in_booz file.", default = False)
	parser.add_argument("--stl", dest="lstl", action='store_true',
		help="Output STL file of VMEC boundary", default = False)
	args = parser.parse_args()
	vmec_wout = VMEC()
	vmec_input = VMEC_INDATA()
	libStell = LIBSTELL()
	if args.vmec_ext:
		linput = False
		loutput = False
		try:
			vmec_input.read_indata('input.'+args.vmec_ext)
			linput = True
		except:
			print(f'Could not file input file: input.{args.vmec_ext}')
		try:
			vmec_wout.read_wout(args.vmec_ext)
			loutput = True
		except:
			print(f'Could not file input file: wout_{args.vmec_ext}.nc or wout.{args.vmec_ext}')
		if not (linput or loutput): sys.exit(-1)
		# Write in_booz file
		if (loutput and args.lbooz):
			filename = 'in_booz.'+args.vmec_ext
			f = open(filename,'w')
			f.write(f'{vmec_wout.mpol*4} {vmec_wout.ntor*2}\n')
			f.write(f'{args.vmec_ext}\n')
			for i in range(vmec_wout.ns):
				f.write(f' {i+1}')
			f.write('\n')
			f.close()
		# Do Input file plot
		if (args.lplot and linput):
			px = 1/pyplot.rcParams['figure.dpi']
			fig=pyplot.figure(figsize=(1024*px,768*px))
			ax=fig.add_subplot(221)
			pyplot.subplots_adjust(hspace=0.4,wspace=0.3)
			s = np.linspace(0.0,1.0,128)
			f = np.zeros(128)
			for i,x in enumerate(s):
				f[i] = libStell.pmass(x)
			ax.plot(s,f/1E3,'k')
			ax.text(0.02,0.19,rf'PHIEDGE={vmec_input.phiedge:4.3f} [Wb]', horizontalalignment='left',\
				verticalalignment='center', transform=ax.transAxes)
			ax.text(0.02,0.12,rf'PRESS_SCALE={vmec_input.pres_scale}', horizontalalignment='left',\
				verticalalignment='center', transform=ax.transAxes)
			ax.text(0.02,0.05,rf'MGRID_FILE: {vmec_input.mgrid_file}', horizontalalignment='left',\
				verticalalignment='center', transform=ax.transAxes)
			ax.text(0.98,0.95,f'PMASS_TYPE: {vmec_input.pmass_type}', horizontalalignment='right',\
				verticalalignment='center', transform=ax.transAxes)
			ax.set_xlim(0,1)
			ax.set_xlabel('Norm. Tor. Flux (s)')
			ax.set_ylabel("Pressure [kPa]")
			ax.set_title(f'VMEC Input: {args.vmec_ext}')
			ax=fig.add_subplot(222)
			if vmec_input.ncurr==1:
				for i,x in enumerate(s):
					f[i] = libStell.pcurr(x)
					ax.set_ylabel("Current I(s)")
					temp_str='PCURR'
					temp_type = vmec_input.pcurr_type
			else:
				for i,x in enumerate(s):
					f[i] = libStell.piota(x)
					ax.set_ylabel(r"$\iota$")
					temp_str='PIOTA'
					temp_type = vmec_input.piota_type
			ax.plot(s,f,'k')
			ax.set_xlabel('Norm. Tor. Flux (s)')
			ax.text(0.98,0.95,f'{temp_str}_TYPE: {temp_type}', horizontalalignment='right',\
				verticalalignment='center', transform=ax.transAxes)
			ax.text(0.02,0.12,rf'CURTOR: {vmec_input.curtor:6.1f}', horizontalalignment='left',\
				verticalalignment='center', transform=ax.transAxes)
			ax.text(0.02,0.05,rf'NCURR: {vmec_input.ncurr}', horizontalalignment='left',\
				verticalalignment='center', transform=ax.transAxes)
			pyplot.show()
		# Do wout file plot
		if (args.lplot and loutput):
			px = 1/pyplot.rcParams['figure.dpi']
			fig=pyplot.figure(figsize=(1024*px,768*px))
			ax=fig.add_subplot(221)
			pyplot.subplots_adjust(hspace=0.4,wspace=0.3)
			ax.plot(np.linspace(0.0,1.0,vmec_wout.ns),vmec_wout.presf/1E3,'k')
			ax.text(0.02,0.47,rf'$B_0$={vmec_wout.b0:4.3f} [T]', horizontalalignment='left',\
				verticalalignment='center', transform=ax.transAxes)
			ax.text(0.02,0.40,rf'$R/a$={vmec_wout.aspect:4.3f}', horizontalalignment='left',\
				verticalalignment='center', transform=ax.transAxes)
			ax.text(0.02,0.33,rf'$R$={vmec_wout.rmajor:4.3f}', horizontalalignment='left',\
				verticalalignment='center', transform=ax.transAxes)
			ax.text(0.02,0.26,rf'$a$={vmec_wout.aminor:4.3f}', horizontalalignment='left',\
				verticalalignment='center', transform=ax.transAxes)
			ax.text(0.02,0.19,rf'$W_p$={vmec_wout.eplasma[0]/1000:4.3f} [kJ]', horizontalalignment='left',\
				verticalalignment='center', transform=ax.transAxes)
			ax.text(0.02,0.12,rf'$<\beta>$={vmec_wout.betatot*100:3.2f} %', horizontalalignment='left',\
				verticalalignment='center', transform=ax.transAxes)
			ax.text(0.02,0.05,rf'MGRID_FILE: {vmec_wout.mgrid_file}', horizontalalignment='left',\
				verticalalignment='center', transform=ax.transAxes)
			ax.text(0.98,0.95,f'PMASS_TYPE: {vmec_wout.pmass_type}', horizontalalignment='right',\
				verticalalignment='center', transform=ax.transAxes)
			ax.set_xlim(0,1)
			ax.set_xlabel('Norm. Tor. Flux (s)')
			ax.set_ylabel("Pressure [kPa]")
			ax.set_title(f'VMEC Ouput: {args.vmec_ext}')
			ax=fig.add_subplot(222)
			ax2 = ax.twinx()
			ax.plot(np.linspace(0.0,1.0,vmec_wout.ns),vmec_wout.iotaf,'k')
			ax2.plot(np.linspace(0.0,1.0,vmec_wout.ns),vmec_wout.jcurv/1E3,'r')
			ax.text(0.05,0.05,rf'$I$={vmec_wout.itor/1000:4.2f} [kA]', horizontalalignment='left',\
				verticalalignment='center', transform=ax.transAxes)
			ax.text(0.98,0.95,f'PCURR_TYPE: {vmec_wout.pcurr_type}', horizontalalignment='right',\
				verticalalignment='center', transform=ax.transAxes)
			ax.text(0.98,0.88,f'PIOTA_TYPE: {vmec_wout.piota_type}', horizontalalignment='right',\
				verticalalignment='center', transform=ax.transAxes)
			ax.set_xlim(0,1)
			ax.set_xlabel('Norm. Tor. Flux (s)')
			ax.set_ylabel(r"$\iota$",color='k')
			ax2.set_ylabel(r"$<j_v>$ [$kA/m^2$]",color='r')
			ax=fig.add_subplot(223)
			theta = np.ndarray((360,1))
			zeta  = np.ndarray((3,1))
			for j in range(360): theta[j]=2.0*np.pi*j/359.0
			for j in range(3):   zeta[j]=     np.pi*j/2.0
			r = vmec_wout.cfunct(theta,zeta,vmec_wout.rmnc,vmec_wout.xm,vmec_wout.xn/vmec_wout.nfp)
			z = vmec_wout.sfunct(theta,zeta,vmec_wout.zmns,vmec_wout.xm,vmec_wout.xn/vmec_wout.nfp)
			ax.plot(r[1,1,0],z[1,1,0],'+r')
			ax.plot(r[1,1,1],z[1,1,1],'+g')
			ax.plot(r[1,1,2],z[1,1,2],'+b')
			j = vmec_wout.ns-1
			ax.plot(r[j,:,0],z[j,:,0],'r')
			ax.plot(r[j,:,1],z[j,:,1],'g')
			ax.plot(r[j,:,2],z[j,:,2],'b')
			j = int(vmec_wout.ns/4)
			ax.plot(r[j,:,0],z[j,:,0],'--r')
			ax.plot(r[j,:,1],z[j,:,1],'--g')
			ax.plot(r[j,:,2],z[j,:,2],'--b')
			ax.set_aspect('equal', adjustable='box')
			ax.text(0.02,0.05,rf'NFP: {vmec_wout.nfp}', horizontalalignment='left',\
				verticalalignment='center', transform=ax.transAxes)
			ax=fig.add_subplot(224)
			theta = np.ndarray((256,1))
			zeta  = np.ndarray((256,1))
			for j in range(256): theta[j]=2.0*np.pi*j/255.0
			for j in range(256):  zeta[j]=2.0*np.pi*j/255.0
			b = vmec_wout.cfunct(theta,zeta,vmec_wout.bmnc,vmec_wout.xm_nyq,vmec_wout.xn_nyq/vmec_wout.nfp)
			j = int(vmec_wout.ns/4)
			h=ax.pcolormesh(np.squeeze(b[j,:,:]),cmap='jet',shading='gouraud')
			ax.set_xlabel(r"$\zeta [rad]$")
			ax.set_ylabel(r"$\theta_{VMEC}$ [rad]")
			ax.set_title("|B| at mid radius")
			fig.colorbar(h,label='[T]')
			pyplot.show()
		# Output an STL file
		if (loutput and args.lstl):
			theta = np.linspace([0],[np.pi*2],512)
			phi   = np.linspace([0],[np.pi*2],512)
			r = vmec_wout.cfunct(theta,phi,vmec_wout.rmnc,vmec_wout.xm,vmec_wout.xn)
			z = vmec_wout.sfunct(theta,phi,vmec_wout.zmns,vmec_wout.xm,vmec_wout.xn)
			vmec_wout.surfaceSTL(r,z,phi,filename='plasma_'+args.vmec_ext+'.stl')
	sys.exit(0)