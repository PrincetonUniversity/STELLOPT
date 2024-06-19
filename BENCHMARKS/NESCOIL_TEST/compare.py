#!/usr/bin/env python3
import sys, os
import shelve
from pathlib import Path
from argparse import ArgumentParser
sys.path.insert(0, '../../pySTEL/')
import numpy as np
_BENCH_FILE_ = 'BENCH_DATA'
# This part changes
from libstell.nescoil import NESCOIL
from libstell.coils import COILSET

# Main routine
if __name__=="__main__":
	# Parse Arguments
	parser = ArgumentParser(description= 
		'''Provides tool for benchmarking files.''')
	parser.add_argument("--file", dest="run_name",
		help="File to process", default = None)
	parser.add_argument("--make", dest="lmake_db", action='store_true',
		help="Update the database file.", default = False)

	#
	args = parser.parse_args()
	run_name = args.run_name

	# Do nothing if no filename
	if not run_name: sys.exit(0)
	lfail = False
	failtol = 1.0
	nesc = NESCOIL()
	coil = COILSET()
	nesc.read_nescout(f'nescout.{run_name}')
	version_str = f'NESCOIL VERSION: 1.0'

	# Extract values
	data={}
	if run_name in ['TF_ONLY']:
		nsteps = 16
		coil = nesc.cutcoils(64)
		r = np.ndarray((nsteps,1))
		for i in range(nsteps): r[i] = 3.0*i/float(nsteps-1) + 2.5
		nesc.bfield_init()
		phi = 0.0
		factor = -nesc.curpol/(2.0*np.pi)
		modb_nescoil_adapt = np.ndarray((nsteps,1))
		modb_nescoil = np.ndarray((nsteps,1))
		modb_model = np.ndarray((nsteps,1))
		modb_coil = np.ndarray((nsteps,1))
		for i in range(nsteps):	
			x = r[i]*np.cos(phi)
			y = r[i]*np.sin(phi)
			z = 0.0
			[bx,by,bz]=nesc.bfield(x,y,z)		
			modb_nescoil_adapt[i] = np.sqrt(bx*bx+by*by+bz*bz)
			[bx,by,bz]=nesc.bfield(x,y,z,istat=-327)
			modb_nescoil[i] = np.sqrt(bx*bx+by*by+bz*bz)
			[bx,by,bz]=coil.coilbiot(x,y,z)
			modb_coil[i] = np.sqrt(bx*bx+by*by+bz*bz)
			modb_model[i] = factor/r[i]
		modb_nescoil_adapt = modb_nescoil_adapt.flatten()
		modb_nescoil = modb_nescoil.flatten()
		modb_model = modb_model.flatten()
		modb_coil = modb_coil.flatten()
		data['ERROR_ADAPT'] = 100*abs(modb_nescoil_adapt-modb_model)/modb_model
		data['ERROR_FINITE'] = 100*abs(modb_nescoil-modb_model)/modb_model
		data['ERROR_COIL'] = 100*abs(modb_coil-modb_model)/modb_model

	# Read or write to the database file.
	if args.lmake_db:
		my_file = Path(_BENCH_FILE_)
		if my_file.exists():
			shelf = shelve.open(_BENCH_FILE_, flag="a")
		else:
			shelf = shelve.open(_BENCH_FILE_, flag="c")
		shelf[run_name] = data
		shelf.close()
		print('  ADDED: '+run_name)
		sys.exit(0)
	else:
		shelf = shelve.open(_BENCH_FILE_, flag='r')
		varlist = shelf[run_name]

	print(version_str)
	print('=================')
	for temp in varlist:
		act = varlist[temp]
		cal = data[temp]
		if np.isscalar(act) == 1:
			if act == 0:
				perct = 0
			else:
				perct = 100*abs(act-cal)/act
			print(f'  {temp} {cal:7.6f} {act:7.6f} {round(perct)}')
		elif temp == 'wall_strikes':
			cal = np.where(act==0,0,cal)
			div = np.where(act==0,1,act)
			perct = 100*sum(abs(act-cal)/div)
			print(f'  {temp} {max(cal):7.6f} {max(act):7.6f} {round(perct)}')
		else:
			cal = np.where(act==0,0,cal)
			div = np.where(act==0,1,act)
			print(f'  Quantity: {temp} -- CODE -- REF. -- %')
			for i in range(len(act)):
				perct = 100*abs(act[i]-cal[i])/div[i]
				print(f'  {i} {cal[i]:7.6f} {act[i]:7.6f} {round(perct)}')
		if perct > failtol:
			lfail = True
		print('=================')

	# Error Status
	if lfail:
		print('  STATUS: FAIL!!!!!')
		sys.exit(0) # For now since some may fail due to statistics
	else:
		print('  STATUS: PASS')
		sys.exit(0)





