# -*- coding: utf-8 -*-

# a small python lib to test reading the jxbout file of a VMEC calculation

from netCDF4 import Dataset
import numpy as np

mu0=4.0e-7*np.pi

rootgrp = Dataset(filename, 'r')
vmec_data=dict()

vmec_data['ns']=rootgrp['/ns'][0]





