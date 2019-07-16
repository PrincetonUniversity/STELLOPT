# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 21:36:05 2019

@author: gawe
"""

def ftell(fid):
    """
        A python function that mimics ftell from MATLAB
            returns the current file position
    """
    return fid.tell()

def fseek(fid, offset, origin='currentposition'):
    """
        A python function that mimics fseek from MATLAB
            sets the current file position
    inputs:
        fid    - file object identifier
        offset - integer number of bytes to move the current file position
        origin - reference from where to offset the file position

    In text files only seeks relative to the beginning of the file are allowed
    (except for seeking to the very end of the file with seek(0, 2))

    The only valid offset values are those returned from fid.tell or 0
    """
    if type(origin) == type(''):
        if origin.lower().find('currentposition')>-1:   # default
            origin = 1
        elif origin.lower().find('start-of-file')>-1:
            origin = 0
        elif origin.lower().find('end-of-file')>-1:
            origin = 2
        # end if
    # end if
    return fid.seek(offset, origin)

def frewind(fid):
    return fseek(fid, 0, 'start-of-file')

def fgets(fid):
    """
        A python function that mimics fgets from MATLAB
            returns the next line of a file while including the new line character
    """
    return fid.readline()

def fgetl(fid):
    """
        A python function that mimics fgetl from MATLAB
            returns the next line of a file while stripping the new line character
    """
    return fgets(fid).replace('\n','')

# ===================================================================== #

def findstr(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub) # use start += 1 to find overlapping matches

# ===================================================================== #
# ===================================================================== #
