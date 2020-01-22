#!/usr/bin/env python3
import sys, os, getopt
import argparse
import numpy as np                    #For Arrays
from libstell.libstell import read_vmec, cfunct, sfunct, torocont, isotoro, calc_jll

try:
	qtCreatorPath=os.environ["STELLOPT_PATH"]
except KeyError:
	print("Please set environment variable STELLOPT_PATH")
	sys.exit(1)

if __name__ == "__main__":
	#print(sys.argv)
	lvmec=0
	lspec=0
	lfocus=0
	lflip=0
	parser = argparse.ArgumentParser(description='Outputs VMEC harmonics to stdout.')
	#print(parser)
	parser.add_argument('filename',help='Filename (wout file) to read.')
	parser.add_argument('-v','--vmec', action='store_true', dest='lvmec', help='Print VMEC harmonics.')
	parser.add_argument('-f','--focus', action='store_true', dest='lfocus', help='Print FOCUS harmonics.')
	parser.add_argument('-s','--spec', action='store_true', dest='lspec', help='Print SPEC harmonics.')
	parser.add_argument('--flip', action='store_true', dest='lflip', help='Flip the jacobian.')
	parser.add_argument('-k', action='store', dest='k', type=int, help='Evaluate Surface (default: ns)')
	args   = parser.parse_args()
	lvmec  = args.lvmec
	lspec  = args.lspec
	lfocus = args.lfocus
	lflip  = args.lflip
	vmec_data=read_vmec(args.filename)
	# note that we use nu+nv in pystel
	k = vmec_data['ns']
	if args.k:
		k = args.k
	# Flip Jacobian
	if lflip:
		for i in range(vmec_data['mnmax']):
			m = vmec_data['xm'][i]
			if m > 0:
				if np.remainder(m,2) == 0:
					vmec_data['zmns'][k-1,i] = -vmec_data['zmns'][k-1,i]
				else:
					vmec_data['rmnc'][k-1,i] = -vmec_data['rmnc'][k-1,i]
		vmec_data['xn'] = - vmec_data['xn']
		print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
		print('!!!!!!!!!!!!!!!! JACOBIAN SIGN FLIPPPED !!!!!!!!!!!!!!!')
		print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
	if lvmec:
		ntor = vmec_data['ntor']
		temp = '  RAXIS_CC = '
		for i in range(ntor+1):
			temp = temp + "{:20.10e}".format(vmec_data['rmnc'][0,i]) + ' '
		print(temp)
		temp = '  ZAXIS_CS = '
		for i in range(ntor+1):
			temp = temp + "{:20.10e}".format(vmec_data['zmns'][0,i]) + ' '
		print(temp)
		if vmec_data['lasym']:
			temp = '  RAXIS_CS = '
			for i in range(ntor+1):
				temp = temp + "{:20.10e}".format(vmec_data['rmns'][0,i]) + ' '
			print(temp)
			temp = '  ZAXIS_CC = '
			for i in range(ntor+1):
				temp = temp + "{:20.10e}".format(vmec_data['zmnc'][0,i]) + ' '
			print(temp)
		for i in range(vmec_data['mnmax']):
			temp = '  RBC(' + "{:3d}".format(int(vmec_data['xm'][i])) + ',' \
			                + "{:3d}".format(-int(vmec_data['xn'][i]/vmec_data['nfp'])) + \
			            ') = ' + "{:20.10e}".format(vmec_data['rmnc'][k-1,i]) + \
			       '    ZBS(' + "{:3d}".format(int(vmec_data['xm'][i])) + ',' \
			                + "{:3d}".format(-int(vmec_data['xn'][i]/vmec_data['nfp'])) + \
			            ') = ' + "{:20.10e}".format(vmec_data['zmns'][k-1,i])
			print(temp)
		if vmec_data['lasym']:	
			for i in range(vmec_data['mnmax']):
				temp = '    RBS(' + "{:3d}".format(int(vmec_data['xm'][i])) + ',' \
				                + "{:3d}".format(-int(vmec_data['xn'][i]/vmec_data['nfp'])) + \
				            ') = ' + "{:20.10e}".format(vmec_data['rmns'][k-1,i]) + \
				       '    ZBC(' + "{:3d}".format(int(vmec_data['xm'][i])) + ',' \
				                + "{:3d}".format(-int(vmec_data['xn'][i]/vmec_data['nfp'])) + \
				            ') = ' + "{:20.10e}".format(vmec_data['zmnc'][k-1,i])
				print(temp)
	if lspec:
		if (vmec_data['lasym']):
			print('mn   m   n   rmnc   zmns   rmns   zmnc')
			for i in range(vmec_data['mnmax']):
				temp = "{:3d}".format(i) + '   ' + \
				       "{:3d}".format(int(vmec_data['xm'][i])) + '   ' + \
				       "{:3d}".format(int(vmec_data['xn'][i]/vmec_data['nfp'])) + '   ' + \
				       "{:20.10e}".format(vmec_data['rmnc'][k-1,i]) + '   ' + \
				       "{:20.10e}".format(vmec_data['zmns'][k-1,i]) + '   ' + \
				       "{:20.10e}".format(vmec_data['rmns'][k-1,i]) + '   ' + \
				       "{:20.10e}".format(vmec_data['zmnc'][k-1,i])
				print(temp)
		else:
			print('mn   m   n   rmnc   zmns   rmns   zmnc')
			for i in range(vmec_data['mnmax']):
				temp = "{:3d}".format(i) + '   ' + \
				       "{:3d}".format(int(vmec_data['xm'][i])) + '   ' + \
				       "{:3d}".format(int(vmec_data['xn'][i]/vmec_data['nfp'])) + '   ' + \
				       "{:20.10e}".format(vmec_data['rmnc'][k-1,i]) + '   ' + \
				       "{:20.10e}".format(vmec_data['zmns'][k-1,i]) + '   ' + \
				       "{:20.10e}".format(0.0) + '   ' + \
				       "{:20.10e}".format(0.0) 
				print(temp)
	if lfocus:
		print('#bmn  bNfp nbf')
		print("{:3d}".format(vmec_data['mnmax']) + ' ' + "{:3d}".format(vmec_data['nfp']) + ' ' + "{:3d}".format(vmec_data['mnmax']))
		print('#------plasma boundary harmonics-------')
		print('# n m Rbc Rbs Zbc Zbs')
		if (vmec_data['lasym']):
			for i in range(vmec_data['mnmax']):
				temp = "{:3d}".format(int(vmec_data['xn'][i]/vmec_data['nfp'])) + '   ' + \
				       "{:3d}".format(int(vmec_data['xm'][i])) + '   ' + \
				       "{:20.10e}".format(vmec_data['rmnc'][k-1,i]) + '   ' + \
				       "{:20.10e}".format(vmec_data['rmns'][k-1,i]) + '   ' + \
				       "{:20.10e}".format(vmec_data['zmnc'][k-1,i]) + '   ' + \
				       "{:20.10e}".format(vmec_data['zmns'][k-1,i])
				print(temp)
		else:
			for i in range(vmec_data['mnmax']):
				temp = "{:3d}".format(int(vmec_data['xn'][i]/vmec_data['nfp'])) + '   ' + \
				       "{:3d}".format(int(vmec_data['xm'][i])) + '   ' + \
				       "{:20.10e}".format(vmec_data['rmnc'][k-1,i]) + '   ' + \
				       "{:20.10e}".format(0.0) + '   ' + \
				       "{:20.10e}".format(0.0) + '   ' + \
				       "{:20.10e}".format(vmec_data['zmns'][k-1,i])
				print(temp)

