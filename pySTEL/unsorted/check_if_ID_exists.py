#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 17:00:20 2017

@author: jons
"""

from osa import Client

vmec = Client("http://esb.ipp-hgw.mpg.de:8280/services/vmec_v6?wsdl")

test_id="aug_asy.example01.01.01"

print(vmec.service.vmecIdentifierExists(test_id))