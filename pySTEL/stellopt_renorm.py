#!/usr/bin/env python3
import sys, os, argparse
import numpy as np                    #For Arrays
from math import pi
from libstell.stellopt import read_stellopt

try:
	qtCreatorPath=os.environ["STELLOPT_PATH"]
except KeyError:
	print("Please set environment variable STELLOPT_PATH")
	sys.exit(1)

if __name__ == "__main__":
	lvmec=0
	lspec=0
	lfocus=0
	lflip=0
	parser = argparse.ArgumentParser(description='Calculates renormaliztion of STELLOPT SIGMAS')
	#print(parser)
	parser.add_argument('filename',help='Filename (stellopt.ext) to read.')
	#parser.add_argument('-v','--vmec', action='store_true', dest='lvmec', help='Print VMEC harmonics.')
	#parser.add_argument('-f','--focus', action='store_true', dest='lfocus', help='Print FOCUS harmonics.')
	#parser.add_argument('-s','--spec', action='store_true', dest='lspec', help='Print SPEC harmonics.')
	#parser.add_argument('--flip', action='store_true', dest='lflip', help='Flip the jacobian.')
	#parser.add_argument('-k', action='store', dest='k', type=int, help='Evaluate Surface (default: ns)')
	#parser.add_argument('-t','--thresh', action='store', dest='thres', type=float, help='Threshold harmonics')
	args   = parser.parse_args()
	failtol = 1.0
	filename='stellopt.BASIC'
	data=read_stellopt(args.filename)
	for name in ['TE','TI','NELINE','TELINE','TILINE','FARADAY','SXR','XICS','XICS_BRIGHT','XICS_W3','XICS_V',
		'B_PROBES','FLUXLOOPS','SEGROG','MSE','ECEREFLECT']:
		if name+'_equil' in data:
			if (np.ndim(data[name+'_equil'])) == 2:
				vals = data[name+'_equil'][0,:]
				sig  = data[name+'_sigma'][0,:]
				targ = data[name+'_target'][0,:]
				chi  = ((targ-vals)/sig)**2
				factor = 1.0/sum(chi)
			else:
				vals = data[name+'_equil'][0]
				sig  = data[name+'_sigma'][0]
				targ = data[name+'_target'][0]
				chi  = ((targ-vals)/sig)**2
				factor = 1.0/chi
			print('!  '+name+': '+str(factor))

sys.exit(0)
