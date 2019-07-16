# -*- coding: utf-8 -*-
"""
  Text and string utilities
"""
# ======================================================================== #
# ======================================================================== #

from __future__ import absolute_import, with_statement, absolute_import, \
                       division, print_function, unicode_literals

__all__ = ['scanf', 'sscanf', 'fscanf']

import numpy as _np

__metaclass__ = type



# ======================================================================== #
# ======================================================================== #
"""
Small MATLAB-esque sscanf-implementation.
    also implemented are ftell, fseek, frewind, fgets, and fgetl

    This implementation does not use regular expressions!
    TODO:!   Upgrade fscanf implementation to mimic the scanf implementation below
"""


nptypes = {'i':_np.int, 'd':_np.int, 'ld':_np.int, 'li':_np.int,  # TODO:! actually use correct datatypes
           'u':_np.uint, 'o':_np.uint, 'x':_np.uint, 'lu':_np.uint, 'lo':_np.uint, 'lx':_np.uint,
           'f':_np.float, 'e':_np.float, 'g':_np.float,
           's':_np.str, 'c':_np.str,
           'b':_np.bool}

#nptypes = {'i':'int', 'd':'int', 'ld':'int', 'li':'int',  # TODO:! actually use correct datatypes
#           'u':'uint', 'o':'uint', 'x':'uint', 'lu':'uint', 'lo':'uint', 'lx':'uint',
#           'f':'float', 'e':'float', 'g':'float',
#           's':'str', 'c':'str',
#           'bool':'bool'}



def scanf(fmt):
    pass

def sscanf(s, fmt):
    pass

#def sscanf(s, fmt, sizeA=None):
#
#    if sizeA is None:
#        return sscanf(s, fmt)
#    # end if
#
#    data = []
#    nrow = 1
#    ncol = 1
#    if len(_np.atleast_1d(sizeA))>1:
#        nrow = sizeA[0]
#        ncol = sizeA[1]
#    else:
#        ncol = sizeA
#    # end if
#    return data

def fscanf(fid, dtypes='%g', size=None, offset=None, eol='\n'):
    """
    This function tries to mimic fscanf from MATLAB. In MATLAB, fscanf
    reads in formatted data from an open file up to a specified shape.
        fscanf(fileID,formatSpec,sizeA)

    Input:
        fid    - file id generated from fid = open(filename), etc.
        dtypes - data types mimicing the formatSpec keyword in MATLAB's fscanf function
            See the note below
        size   - 1D or 2D shape of returned data, [nrow, ncol] or an integer
            The shape of the output data
            if size = 13, then 13 elements will be read in along a single row
            if size = [3, 2], then 2 elements will be read in for three rows
        offset - The starting point to read from the file
            MATLAB automatically keeps track.
        unpack - output a tuple for expansion into multiple variables
            default to False, automatically True if multiple data types are entered
    Outputs:
        data - shape is equal to input size
                if multiple data types are input
                    output structure is a list
                if single data type requested and more than 1 element:
                    output structure is a numpy array
                caveat:  if EOF is encountered during reading a shaped array,
                        then a list of lists will be output

    MATLAB cheat-sheet for fscanf

    Formatting
        Integer, signed
            %d - Base 10
            %i - The values in the file determine the base:
                    The default is base 10.
                    If the initial digits are 0x or 0X, then the values are hexadecimal (base 16).
                    If the initial digit is 0, then values are octal (base 8).
            %ld or %li - 64-bit values, base 10, 8, or 16
        Integer, unsigned
            %u - Base 10
            %o - Base 8 (octal)
            %x - Base 16 (hexadecimal)
            %lu, %lo, %lx - 64-bit values, base 10, 8, or 16
        Floating-point number
            %f, %e, %g - Floating-point fields can contain any of the following
                (not case sensitive): Inf, -Inf, NaN, or -NaN.

        String
            %s - string - Read all characters excluding white spaces
            %c - string - Read any single character, including white space.
                    To read multiple characters at a time, specify field width.
    """
    # ===================== #
    if offset is not None:
        fid.seek(offset, 0)
    # end if
#    dtypes_in = dtypes

    # ===================== #
    # Parse the input data types for the delimiter
    # Convert the different datatypes into a list
    delim = ' '
    if dtypes.find(',')>-1:   # delimiter hidden in fmt string
        delim = ','
        dtypes = dtypes.replace(',', '')
    elif dtypes.find(';')>-1:   # delimiter hidden in fmt string
        delim = ';'
        dtypes = dtypes.replace(';', '')
    elif dtypes.find('\t')>-1:   # delimiter hidden in fmt string
        delim = '\t'
        dtypes = dtypes.replace('\t', '')
    # end if

    kwargs = {'delim':delim, 'eol':eol}
    if dtypes == '\n':
        # special request to read to end of line and stop there
        return eztxt.read_rol(fid, **kwargs)

    # ===================== #
    # Convert to python / numpy data types
    global nptypes

    # Count and parse the requested data types input by the user
    dtypes = list(_np.atleast_1d(dtypes.split('%')))[1:]
    ndt = len(dtypes)
    dtypes = [nptypes[dtypes[ii]] for ii in range(ndt)]

    # ===================== #
    # now parse the size keyword input for easier handling later
    nrow = 1
    ncol = 1
    if size is not None:
        size = _np.atleast_1d(size)
        if len(size)>1:
            nrow = size[0]
            ncol = size[1]
        else:
            ncol = size[0]
        # end if
    # end if

    # ===================== #

    if ndt>1 and (dtypes[0] != _np.atleast_1d(dtypes)).any():
        # multiple datatypes entered, we have to assign data types after getting data
        # or perform multiple calls
        kwargs['dtype'] = _np.str
    else:
        # single data type entered and numpy/column vectors requested
        # we can allow the reader to do the data typing on the fly
        kwargs['dtype'] = dtypes[0]
    # end if

    if nrow>1:
        # if an array was input for size then read in to the specified shape
        reader = eztxt.read_to_shape
        args = (fid, size)
    else:
        if ncol>1:
            # if more than one element was requested
            reader = eztxt.read_to_count
            args = (fid, ncol)
        else:
            # if a single-element was requested
            reader = eztxt.read_element
            args = (fid,)
        # end if
    # end if

    # Get the data and return it as a list (or list of lists)
    try:
        data = reader(*args, **kwargs)
    except EOFError:
        print('premature termination by end-of-file (before input size was reached)')
    except:
        raise
    # end try

    # get actual output dimensions returned
    sh = _np.shape(_np.atleast_1d(data))
    if len(sh)>1:
        nr, nc = sh
    else:
        nr = 1
        nc = sh[0]
    # end if

    # ======================= #
    if ndt>1 and (dtypes[0] != _np.atleast_1d(dtypes)).any():
        # Loop over the elements of the list and assign the data type
        if nr > 1: # nrow
            for ii in range(nr):  # nrow
                nc = len(_np.atleast_1d(data[ii]))
                for jj in range(nc):   # ncol
                    data[ii][jj] = dtypes[ii](data[ii][jj])
                # end for
            # end for
        else:
            data = [dtypes[jj](data[jj]) for jj in range(nc)] # ncol
        # end if
    else:
        if nr == 1:
            if dtypes[0] == str:
                # for reading in strings that have spaces in them (when delim==' ')
                if data is None:
                    data = ' '
                data = ''.join(data)
                data = data.strip(delim)
            elif nc > 1:
                data = _np.asarray(data)
            # end if
        elif nr>1:
            data = _np.asarray(data)
        # end if
    # end if
    # ======================= #
    return data

#def fscanf(fid, dtypes='%g', size=None, offset=None, unpack=False, eol='\n'):
#    """
#    This function tries to mimic fscanf from MATLAB. In MATLAB, fscanf
#    reads in formatted data from an open file up to a specified shape.
#        fscanf(fileID,formatSpec,sizeA)
#
#    Input:
#        fid    - file id generated from fid = open(filename), etc.
#        dtypes - data types mimicing the formatSpec keyword in MATLAB's fscanf function
#            See the note below
#        size   - 1D or 2D shape of returned data, [nrow, ncol] or an integer
#            The shape of the output data
#            if size = 13, then 13 elements will be read in along a single row
#            if size = [3, 2], then 2 elements will be read in for three rows
#        offset - The starting point to read from the file
#            MATLAB automatically keeps track.
#        unpack - output a tuple for expansion into multiple variables
#            default to False, automatically True if multiple data types are entered
#    Outputs:
#        data - shape is equal to input size, but defaults to (1,length of the row)
#
#    """
#    kwargs = {'delim':delim, 'eol':eol}
#    if ndt>1:
#        # multiple datatypes entered, we have to assign data types after getting data
#        # or perform multiple calls
#        kwargs['dtype'] = _np.str
#
#        # if multiple data formats have been entered we have to unpack it
#        unpack = True
#    else:
#        # single data type entered and numpy/column vectors requested
#        # we can allow the reader to do the data typing on the fly
#        kwargs['dtype'] = dtypes[0]
#    # end if
#
#    if nrow>1:
#        # if an array was input for size then read in to the specified shape
#        reader = eztxt.read_to_shape
#        args = (fid, size)
#    else:
#        if ncol>1:
#            # if more than one element was requested
#            reader = eztxt.read_to_count
#            args = (fid, size)
#        else:
#            # if a single-element was requested
#            reader = eztxt.read_element
#            args = (fid,)
#        # end if
#    # end if
#    try:
#        data = reader(*args, **kwargs)
#    except EOFError:
#        print('premature termination by end-of-file (before input size was reached)')
#    except:
#        raise
#    # end try
#
#    # ======================= #
#    if ndt>1:
#        # Loop over the elements of the list and assign the data type
#        if nrow > 1:
#            for ii in range(nrow):
#                for jj in range(ncol):
#                    data[ii][jj] = dtypes[jj](data[ii][jj])
#                # end for
#            # end for
#        else:
#            for jj in range(ncol):
#                data[jj] = dtypes[jj](data[jj])
#            # end for
#        # end if
#    # end if
#    # ======================= #
#
#    if nrow > 1 and ndt == 1:
#        data = _np.asarray(data)
#        if unpack:
#            # single data format input, column vector output required
#            #  data = fscanf(fid, dtype='%g', [3, 5], unpack=True)
#            #  returns 5 column vectors of float64 data, size(3,)
#            return tuple([data[:,ii] for ii in range(ncol)])
#            # return tuple(map(None, zip(*data)))
#        else:
#            # single data format input, array output required
#            #  data = fscanf(fid, dtype='%g', [3, 5])
#            #  Should output a numpy array sized 3,5 of float64 data
#            return data
#        # end if
#    elif nrow >1 and ndt > 1:
#        # We automatically have to unpack the data in this case
#        #   (unpack == True when ndt>1)
#        # multiple data formats input, expecting column vectors output
#        #   a, b, c = fscanf(fid, dtype='%d%d%g', [100, 3])
#        #Should return an unzipped tuple for expansion
#        #   size(a) == size(b) == size(c) == [100,]
#        return tuple(map(_np.asarray, zip(*data)))
#    elif nrow == 1 and ndt == 1:
#        if unpack:
#            # single data format input, expanded output required
#            #  a, b, c = fscanf(fid, dtype='%i', 3)
#            # Should return a tuple for expansion
#            if dtypes[0] == str:
#                data = data.strip(delim)
#            # end if
#            return tuple(data)
#        else:
#            if dtypes[0] == str:
#                # for reading in strings that have spaces in them (when delim==' ')
#                data = ''.join(data)
#                data = data.strip(delim)
#            # end if
#            return data
#    elif nrow == 1 and ndt > 1:
#       # We automatically have to unpack the data in this case
#        #   (unpack == True when ndt>1)
#        return tuple(data)
#
#        #    multiple data formats input, unzipped tuple output required
#        #        a, b, c = fscanf(fid, dtype='%d%d%g', 3)
#        #            Should return an unzipped tuple for expansion
#        return tuple([nptypes[dtypes[jj]](data[0][jj]) for jj in range(ncol)])
#
#        elif ndt == ncol and nrow>1:
#            #    multiple data formats input, expecting column vectors output
#            #        a, b, c = fscanf(fid, dtype='%d%d%g', [100, 3])
#            #            Should return an unzipped tuple for expansion
#            #        size(a) == size(b) == size(c) == [100,]
#            # Loop over the elements of the list and assign the data type
#            for ii in range(nrow):
#                for jj in range(ncol):
#                    data[ii][jj] = nptypes[dtypes[jj]](data[ii][jj])
#                # end for
#            # end for
#            return tuple(map(_np.asarray, zip(*data)))
#        # end if
#    # end if
#
#    # we will get to this situation when it arises
#    print('unknown situation ... outputing a list of lists from fscanf')
#    print('Note:  no data formatting applied: string output!')
#    return data
#    # end if
#    return data, fid.tell()


class eztxt(object):
    def __init__(self, fid, dtype=str, delim=' ', eol='\n'):
        self.fid = fid           # file id
        self.delim = delim       # delimiter
        self.eol = eol           # line feed
        self.cr                  # carriage return
        self.dtype = dtype       # function handle for data type conversion
    # end def

    @staticmethod
    def read_rol(fid, dtype=str, delim=' ', eol='\n', allowed_white_lines=100):
        buffsize = max((len(delim),len(eol), 1))
        tst = ''
        white_lines = 0
        while True and white_lines<allowed_white_lines:
            tmp = fid.read(buffsize)
            if len(tmp) == 0:
                white_lines += 1  # same result for eof and a blank line
                continue
            # end if
            tst += tmp

            if tst == delim:  # remove new line if it is the only element
                tst = tst.strip(delim)
                continue
            if tst.find(eol)>-1:
                break
            # end if
        # end while
        return tst

#    @staticmethod
#    def read_element(fid, dtype=str, delim=' ', eol='\n', allowed_white_lines=1000, readline=False):
#        tst = _np.fromfile(fid, dtype=dtype, count=1, sep=delim.strip())
#        return tst[0]

    @staticmethod
    def read_element(fid, dtype=str, delim=' ', eol='\n', allowed_white_lines=10, readline=False):
        buffsize = max((len(delim),len(eol), 1))
        tst = ''
        white_lines = 0
        eof = False
#        pos = -1
        while not eof and white_lines<allowed_white_lines:
#            pos = ftell(fid)
            tmp = fid.read(buffsize)

#            if ftell(fid) == pos:
            if len(tmp) == 0:
                white_lines += 1  # same result for eof and a blank line
                continue
            if white_lines>allowed_white_lines:
                # same position before and after read means that we have hit the end-of-file
                eof = True
                break
            # end if

            # ================================== #
            # Append the data
            tst += tmp
#            print(tst)   # debugging makes pretty pyramid shapes on your stdout

            # ================================== #
            if tst == eol:  # remove new line if it is the only element
                tst = tst.strip(eol)
                continue
            if tst == delim:  # remove new line if it is the only element
                tst = tst.strip(delim)
                continue

            if tst.find(delim)>-1 and (not readline):
                if delim != ' ':
                    tst = tst.strip(delim)
                break
            if tst.find(eol)>-1:
                tst = tst.strip(eol)
                break
            # end if
            # ================================== #
        # end while
        if len(tst) == 0:
            return None
        else:
            return dtype(tst)
        # end if

    @staticmethod
    def read_line(fid, dtype=str, delim=' ', eol='\n', allowed_white_lines=100):
        return eztxt.read_element(fid, dtype=dtype, delim=delim, eol=eol, readline=True)

    @staticmethod
    def read_to_count(fid, count, dtype=str, delim=' ', eol='\n'):
        """
        returns a list with length count
        """
        tst = []
        for ii in range(count):
            tmp = eztxt.read_element(fid, dtype=dtype, delim=delim, eol=eol)
            if tmp is not None:
                tst.append(tmp)
            # end if
        # end for
        return tst
#        tst = _np.fromfile(fid, dtype=dtype, count=count, sep=delim.strip())
#        return tst.tolist()

    @staticmethod
    def read_to_shape(fid, shape, dtype=str, delim=' ', eol='\n'):
        """
        returns a list of lists with shape determined by input shape
        """
        nrow = 1
        if len(shape)>1:
            nrow = shape[0]
            ncol = shape[1]
        else:
            ncol = shape[0]
        # end if
        nrow = int(nrow)
        ncol = int(ncol)

        tst = []
        for ii in range(nrow):
            tmp = eztxt.read_to_count(fid, count=ncol, dtype=dtype, delim=delim, eol=eol)
            if (tmp is not None):
                tst.append(tmp)
            # end if
        # end for
        return tst
    # end def methods
# end class eztxt


# ===================================================================== #
# ===================================================================== #

def __test_simple_save():
    tstfil = 'utils.testfile'
    npbool = True
    npfloatarray = _np.ones((3,4), dtype=float)
    npfloatarray += _np.random.normal(0, 1, _np.shape(npfloatarray))
    npint1 = 134
    npint2 = 323
    npint3 = 51
    npfloat = 23.3
    strtst1 = 'what is my name?'
    strtst2 = 'wout.txt'
    tst3 = [3e-4, 1.56e-3]
    testdat = [npbool, npfloatarray, npint1, [npint1, npint2, npint3, npfloat], strtst1, strtst2, tst3]

    eol = '\n'
    delim = ' '
    try:
        fid = open(tstfil, 'w+')
        for tst in testdat:
            if len(_np.atleast_1d(tst))>1:
                try:    tst = tst.flatten()
                except: pass
                # end try

                for ii in range(len(tst)):
                    if ii>0 and ii%4 == 0:  # maximum of 4 elements per line
                        fid.write(eol)
                    # end if
                    fid.write(str(tst[ii])+delim)
                # end for
#                fid.write(' '.join(str(tst[ii]) for ii in range(len(tst))))
            else:
                fid.write(''.join(str(tst)) + delim + eol)
            # end if
            fid.write(eol)  # the write statement just fills a buffer ...
            # end if

            # writing to file only happens when the close statement is called.
            fid.flush()  # forces the buffer to empty into the file
        # end for
    except:
        raise
    finally:
        try:        fid.close()
        except:     pass
    # end try
    return testdat, tstfil
# end def

def test_fscanf():
    import os as _os
    testdat, tstfil = __test_simple_save()

#    elmnts_shape = [1, [3,4], 1, 4, 4, 1, 2]  # actual shapes
#    dtypes = [_np.bool, _np.float64, _np.int, _np.int, _np.str, _np.float64]

    elmnts_shape = [1, [3,4], 1, 4, 4, 1, [1,3]]  # intended to trip the EOFError and handle it appropriately
    dtyps = ['%b','%g','%d','%d%d%d%g','%s','%s', '%e']
    try:
        fid = open('utils.testfile', 'r')

#        kwargs = {'delim':delim, 'eol':eol}
        for ii in range(len(elmnts_shape)):
            sh = elmnts_shape[ii]
            tst = fscanf(fid, dtypes=dtyps[ii], size=sh)

            # test if the data reading worked
            if len(_np.atleast_1d(testdat[ii]))>1 or len(_np.atleast_1d(tst))>1:
                if _np.atleast_1d(tst == testdat[ii]).all():
                    print((tst, 'success'))
                else:
                    print((tst, testdat[ii], 'fail'))
            elif (tst == testdat[ii]):
                print((tst, 'success'))
            else:
                print((tst, testdat[ii], 'fail'))
        # end for
    except:
        raise
    finally:
        fid.close()
    # end try

    _os.remove(tstfil)
# end def

# ===================================================================== #


if __name__ == "__main__":
    import doctest
    test_fscanf()
    doctest.testmod(verbose=True, report=True)
# end if

# ===================================================================== #
# ===================================================================== #

