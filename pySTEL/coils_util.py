#!/usr/bin/env python3
# -*- coding: utf-8 -*-


if __name__=="__main__":
	import sys
	from argparse import ArgumentParser
	from libstell.coils import COILSET
	from libstell.wall import WALL
	from libstell.plot3D import PLOT3D
	from libstell.libstell import FourierRep
	import matplotlib.pyplot as pyplot
	import numpy as np
	from datetime import datetime
	from stl import mesh
	parser = ArgumentParser(description= 
		'''Provides class for accessing coils files also serves as a
		   simple tool for assessing coils or coils files.''')
	parser.add_argument("-c", "--coil", dest="coils_file",
		help="Coils file for input", default = None)
	parser.add_argument("-p", "--plot", dest="lplot", action='store_true',
		help="Plot the coils file.", default = False)
	parser.add_argument("-prz", "--plotRZ", dest="lplotRZ", action='store_true',
		help="Plot each coil group in RZ.", default = False)
	parser.add_argument("--plotcoildist", dest="lplotcoildist", action='store_true',
		help="Plot the coil-coil distance.", default = False)
	parser.add_argument("--plotvolcoil", dest="heightwidth",
		help="Plot volumetric coil of given width and height [m].", default = None)
	parser.add_argument("-b", "--bfield", dest="bxyz",
		help="Output B field at x,y,z", default = None)
	parser.add_argument("-a", "--afield", dest="axyz",
		help="Output A field at x,y,z", default = None)
	parser.add_argument("--genwall", dest="wall_offset",
		help="Generate a wall file based on offset in [m]", default = None)	
	parser.add_argument("--new_pts", dest="new_pts",
		help="Spline coils to new_pts number of points", default = None, type = int)
	parser.add_argument("--fit_surf", dest="lfit_surf", action='store_true',
		help="Fit a surface to a coil.", default = False)
	parser.add_argument("-o", "--output", dest="loutput", action='store_true',
		help="Output the coil", default = False)
	parser.add_argument("--stl", dest="heightwidth_stl",
		help="Generate STL of coil of given width and height [m].", default = None)
	args = parser.parse_args()
	coils = COILSET()
	if args.coils_file: 
		coils.read_coils_file(args.coils_file)
		if args.new_pts: coils.rescalecoils(args.new_pts)
		if args.lplot: coils.plotcoils()
		if args.lplotRZ: coils.plotcoilsRZ()
		if args.lplotcoildist:
			coils.coilCoilDist()
			coils.plotcoilcoilDist()
		if args.heightwidth:
			height,width = args.heightwidth.split(',')
			[vertex,faces] = coils.blenderCoil(height=float(height),width=float(width),lfield_period=True)
			vol_coil = WALL()
			vol_coil.vertex = np.array(vertex)
			vol_coil.faces = np.array(faces, dtype=int)
			vol_coil.nfaces = vol_coil.faces.shape[0]
			vol_coil.nvertex = vol_coil.vertex.shape[0]
			vol_coil.plot_wall_3D()
		if args.heightwidth_stl:
			height,width = args.heightwidth_stl.split(',')
			[vertex,faces] = coils.blenderCoil(height=float(height),width=float(width),lfield_period=False)
			vertex = np.array(vertex)
			faces  = np.array(faces, dtype=int)
			print(vertex.shape)
			print(faces.shape)
			nfaces = faces.shape[0]
			wall_mesh = mesh.Mesh(np.zeros(nfaces, dtype=mesh.Mesh.dtype))
			for i, f in enumerate(faces):
				for j in range(3):
					wall_mesh.vectors[i][j] = vertex[f[j],:]
			wall_mesh.save(args.coils_file+'.stl')
		if args.loutput: coils.write_coils_file(args.coils_file+'_new')
		if args.axyz:
			x,y,z = args.axyz.split(',')
			ax,ay,az = coils.coilvecpot(float(x),float(y),float(z))
			print(f"Vector Potential ({x},{y},{z}) : {ax}, {ay}, {az} ")
		if args.bxyz:
			x,y,z = args.bxyz.split(',')
			bx,by,bz = coils.coilbiot(float(x),float(y),float(z))
			print(f"B-Field ({x},{y},{z}) : {bx}, {by}, {bz} [T]")
		if args.wall_offset:
			vertex = coils.coiloffset(float(args.wall_offset))
			wall = WALL()
			wall.genWallfromOrdered(vertex)
			wall.name = f"Offset from Coils file {args.wall_offset} {args.coils_file}"
			print(wall.vertex.shape)
			wall.write_wall('wall_test.dat')
			# VTK stuff
			renderer = vtk.vtkRenderer()
			render_window = vtk.vtkRenderWindow()
			render_window.AddRenderer(renderer)
			render_window_interactor = vtk.vtkRenderWindowInteractor()
			render_window_interactor.SetRenderWindow(render_window)
			render_window.SetSize(1024, 768)
			wall.plot_wall_3D(renderer=renderer,render_window=render_window)
			coil.plotcoils(renderer=renderer,render_window=render_window)
			render_window.Render()
			render_window_interactor.Start()
		if args.lfit_surf:
			# VTK stuff
			plt3d = PLOT3D()
			print('  Fitting Surface')
			coils.plotcoilsHalfFP(plot3D=plt3d)
			xm,xn,rmnc,zmns=coils.fitSurface()
			theta = np.linspace([0],[np.pi*2],64)
			phi   = np.linspace([0],[np.pi*2],64)/coils.nfp
			FR = FourierRep()
			f=open('coil_surf_harmonics.csv','w')
			mnmax = xn.shape[0]
			f.write(f'# n m Rbc Rbs Zbc Zbs\n')
			for mn in range(mnmax):
				f.write(f'{int(xn[mn,0])}, {int(xm[mn,0])}, {rmnc[0,mn]:10.9e}, {0.0:10.9e}, {0.0:10.9e}, {zmns[0,mn]:10.9e}\n')
			f.close()
			r = FR.cfunct(theta,phi,rmnc,xm,xn)
			z = FR.sfunct(theta,phi,zmns,xm,xn)
			FR.isotoro(r,z,phi,0,plot3D=plt3d,lclosev=False)
			plt3d.render()

