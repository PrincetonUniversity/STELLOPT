#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Main routine
if __name__=="__main__":
	import sys
	from argparse import ArgumentParser
	import numpy as np
	import matplotlib.pyplot as pyplot
	from libstell.vmec import VMEC
	from libstell.beams3d import BEAMS3D_INPUT
	from libstell.plasma import EC,DA
	parser = ArgumentParser(description= 
		'''Utility for generating a BEAMS3D input namelist.''')
	parser.add_argument("-v", "--vmec", dest="vmec_ext",
		help="VMEC file extension", default = None)
	args = parser.parse_args()
	beams3d_input = BEAMS3D_INPUT()
	beams3d_input.read_input('')
	E = 3.5E6*EC
	M = 4.002603*DA
	npitch = 16
	vllov = np.linspace(0.01,0.5,npitch)
	if args.vmec_ext:
		vmec_data = VMEC()
		vmec_data.read_wout(args.vmec_ext)
		beams3d_input.nr = 128
		beams3d_input.nz = 128
		beams3d_input.nphi = int(360/vmec_data.nfp)
		dr = (vmec_data.rmax_surf-vmec_data.rmin_surf)*0.5
		R0 = (vmec_data.rmax_surf+vmec_data.rmin_surf)*0.5
		beams3d_input.rmin = R0 - dr*1.2
		beams3d_input.rmax = R0 + dr*1.2
		beams3d_input.zmin =-vmec_data.zmax_surf*1.2
		beams3d_input.zmax = vmec_data.zmax_surf*1.2
		beams3d_input.phimin = 0.0
		beams3d_input.phimax = 2.0*np.pi/vmec_data.nfp
		beams3d_input.int_type = 'LSODE'
		beams3d_input.follow_tol = 1.0E-8
		beams3d_input.npoinc = 1000
		beams3d_input.vc_adapt_tol = 1.0E-3
		beams3d_input.ns_prof1 = 16
		beams3d_input.ns_prof2 = 2
		beams3d_input.ns_prof3 = 2
		beams3d_input.ns_prof4 = 2
		beams3d_input.ns_prof5 = 2
		beams3d_input.partvmax = float(np.sqrt(E*2.0/M)*1.1)
		# Now determine some things
		sdex = np.linspace(0,vmec_data.ns-1,8,dtype=int)
		theta = np.linspace([0],[np.pi*2],256)
		phi   = np.linspace([0],[np.pi*2],256)/vmec_data.nfp
		r = vmec_data.cfunct(theta,phi,vmec_data.rmnc,vmec_data.xm,vmec_data.xn)
		z = vmec_data.sfunct(theta,phi,vmec_data.zmns,vmec_data.xm,vmec_data.xn)
		b = vmec_data.cfunct(theta,phi,vmec_data.bmnc,vmec_data.xm_nyq,vmec_data.xn_nyq)
		r_start_in = []; z_start_in = []; phi_start_in = []
		vll_start_in = []; mu_start_in = []
		for k in sdex:
			temp=np.argwhere(b[k,:,:] == np.min(b[k,:,:]))
			l = temp[0][0]
			m = temp[0][1]
			r_temp = r[k,l,m]
			z_temp = z[k,l,m]
			b_temp = b[k,l,m]
			p_temp = phi[m]
			v_temp = np.sqrt(2.0*E/M)
			vll_temp = v_temp*vllov
			mu_temp  = 0.5*(v_temp*v_temp-vll_temp*vll_temp)*M/b_temp
			r_start_in.extend([r_temp]*npitch)
			z_start_in.extend([z_temp]*npitch)
			phi_start_in.extend([p_temp]*npitch)
			vll_start_in.extend([vll_temp])
			mu_start_in.extend([mu_temp])
		beams3d_input.r_start_in   = np.array(r_start_in[:]).flatten()
		beams3d_input.z_start_in   = np.array(z_start_in).flatten()
		beams3d_input.phi_start_in = np.array(phi_start_in).flatten()
		beams3d_input.vll_start_in = np.array(vll_start_in).flatten()
		beams3d_input.mu_start_in  = np.array(mu_start_in).flatten()
		beams3d_input.charge_in    = np.ones(len(r_start_in))*EC*2.0
		beams3d_input.mass_in      = np.ones(len(r_start_in))*M
		beams3d_input.zatom_in     = np.ones(len(r_start_in))*2.0
		beams3d_input.t_end_in     = np.ones(len(r_start_in))*100E-3
		beams3d_input.nparticles_start = len(r_start_in)
		beams3d_input.write_input('input.'+args.vmec_ext)



