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
	parser.add_argument("--fullorbit", dest="lfullorbit", action='store_true',
		help="Create Gyro Orbit (full orbit) run.", default = None)
	args = parser.parse_args()
	beams3d_input = BEAMS3D_INPUT()
	beams3d_input.read_input('')
	E = 3.5E6*EC
	M = 4.002603*DA
	npitch = 8
	temp = np.linspace(0.01,0.5,npitch)
	vllov = np.append(-temp[-1::-1],temp)
	npitch = vllov.shape[0]
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
		vr_start_in = []; vphi_start_in = []; vz_start_in = []
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
			vperp_temp = np.sqrt(v_temp*v_temp-vll_temp*vll_temp)
			mu_temp  = 0.5*vperp_temp*vperp_temp*M/b_temp
			if args.lfullorbit:
				for vperp in vperp_temp:
					r_temp = np.squeeze(r[k,l,m])
					z_temp = np.squeeze(z[k,l,m])
					p_temp = np.squeeze(phi[m])
					rg = M*vperp/(2.0*EC*b_temp)
					br,bphi,bz,s,u,info = vmec_data.getBcyl(r_temp,p_temp,z_temp)
					bx = br * np.cos(p_temp) - bphi * np.sin(p_temp)
					by = br * np.sin(p_temp) + bphi * np.cos(p_temp)
					bx = np.squeeze(bx / b_temp)
					by = np.squeeze(by / b_temp)
					bz = np.squeeze(bz / b_temp)
					bpx = -bx * bz * rg
					bpy = -by * bz * rg
					bpz =  (bx * bx + by * by) * rg
					rot_matrix = np.zeros((3,3))
					# theta = 0 solution
					rot_matrix[0,0] = 1.0 
					rot_matrix[0,1] = -bz
					rot_matrix[0,2] =  by
					rot_matrix[1,0] =  bz
					rot_matrix[1,1] = 1.0
					rot_matrix[1,2] = -bx
					rot_matrix[2,0] = -by
					rot_matrix[2,1] =  bx
					rot_matrix[2,2] = 1.0
					x_temp = np.matmul(rot_matrix,[bpx,bpy,bpz])
					xg = by*x_temp[2]-bz*x_temp[1]
					yg = bz*x_temp[0]-bx*x_temp[2]
					zg = bx*x_temp[1]-by*x_temp[0]
					rg = 1.0/(xg*xg+yg*yg+zg*zg)
					x_temp[0] = x_temp[0] + r_temp * np.cos(p_temp)
					x_temp[1] = x_temp[1] + r_temp * np.sin(p_temp)
					x_temp[2] = x_temp[2] + z_temp
					r_temp = np.sqrt(x_temp[0] * x_temp[0] + x_temp[1] * x_temp[1])
					p_temp = np.arctan2(x_temp[1],x_temp[0])
					z_temp = x_temp[2]
					vx_temp = vll_temp*bx + vperp_temp * xg * rg
					vy_temp = vll_temp*by + vperp_temp * yg * rg
					vz_temp = vll_temp*bz + vperp_temp * zg * rg
					vr_temp   = vx_temp*np.cos(p_temp) + vy_temp * np.sin(p_temp)
					vphi_temp =-vx_temp*np.sin(p_temp) + vy_temp * np.cos(p_temp)
					vr_start_in.extend([vr_temp])
					vphi_start_in.extend([vphi_temp])
					vz_start_in.extend([vz_temp])
					r_start_in.extend([r_temp])
					z_start_in.extend([z_temp])
					phi_start_in.extend([p_temp])
			else:
				r_start_in.extend([r_temp]*npitch)
				z_start_in.extend([z_temp]*npitch)
				phi_start_in.extend([p_temp]*npitch)
			vll_start_in.extend([vll_temp])
			mu_start_in.extend([mu_temp])
		beams3d_input.r_start_in   = np.array(r_start_in).flatten()
		beams3d_input.z_start_in   = np.array(z_start_in).flatten()
		beams3d_input.phi_start_in = np.array(phi_start_in).flatten()
		beams3d_input.vll_start_in = np.array(vll_start_in).flatten()
		beams3d_input.mu_start_in  = np.array(mu_start_in).flatten()
		beams3d_input.charge_in    = np.ones(len(r_start_in))*EC*2.0
		beams3d_input.mass_in      = np.ones(len(r_start_in))*M
		beams3d_input.zatom_in     = np.ones(len(r_start_in))*2.0
		beams3d_input.t_end_in     = np.ones(len(r_start_in))*100E-3
		beams3d_input.nparticles_start = len(r_start_in)
		if args.lfullorbit:
			beams3d_input.rho_fullorbit = 0.0
			beams3d_input.vr_start_in   = np.array(vr_start_in[:]).flatten()
			beams3d_input.vphi_start_in   = np.array(vphi_start_in[:]).flatten()
			beams3d_input.vz_start_in   = np.array(vz_start_in[:]).flatten()
		beams3d_input.write_input('input.'+args.vmec_ext)
	sys.exit(0)



