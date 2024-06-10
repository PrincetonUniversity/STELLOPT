#!/usr/bin/env python3
# -*- coding: utf-8 -*-


if __name__=="__main__":
	import sys
	from argparse import ArgumentParser
	from libstell.nescoil import NESCOIL
	import matplotlib.pyplot as pyplot
	import numpy as np
	from datetime import datetime
	parser = ArgumentParser(description= 
		'''Provides class for accessing nescoil files''')
	parser.add_argument("--output", dest="nescout_file",
		help="NESCOIL output file", default = None)
	parser.add_argument("-p", "--plot", dest="lplot", action='store_true',
		help="Plot the coils file.", default = False)
	args = parser.parse_args()
	nescout = NESCOIL()
	if args.nescout_file: 
		nescout.read_nescout(args.nescout_file)
		if args.lplot: nescout.plotpotential()
	sys.exit(0)
