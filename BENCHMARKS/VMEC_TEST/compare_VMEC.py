#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Caoxiang Zhu (czhu@ppp.gov)
For any help, type ./compare_VMEC.py -h
"""
import numpy as np
import xarray
import argparse

# parse command line arguments
# =============================================================================
parser = argparse.ArgumentParser(description="Compare two VMEC netcdf outputs")
parser.add_argument("filename", type=str, help="file name to be compared")
parser.add_argument("reference", type=str, help="reference data")
parser.add_argument("-t", "--tol", type=float, default=1E-12, help="difference tolerance")
 
args = parser.parse_args()
print('Compare VMEC outputs in {:s} and {:s} with tolerance {:12.5E}'.format(args.filename, args.reference, args.tol))
data = xarray.open_dataset(args.filename)
ref = xarray.open_dataset(args.reference)

tol = args.tol
match = True
skip_list = ['version_']

print('======================================')
for key in ref:
    if key in skip_list:
        continue
    try:
        diff = np.linalg.norm(ref[key].values - data[key].values)
        unmatch = diff > tol
        if unmatch:
            match = False
            print('UNMATCHED: '+key, ', diff={:12.5E}'.format(diff))
    except TypeError:
        pass
print('  comparison passed: {:} '.format(match))
print('======================================')
assert match, 'Differences in some elements are larger than the tolerence.'

exit
