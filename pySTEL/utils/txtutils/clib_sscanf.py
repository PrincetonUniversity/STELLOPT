# -*- coding: utf-8 -*-
"""
The code automatically searches for an appropriate library on your path
    on linux:
        ctypes.util.find_library('c')
        should return 'libc.so' or something similar etc.
    on mac:
        ctypes.util.find_library('c')
        should return the full path to libc
    on windows we use the standard distributable
        ctypes.util.find_library('MSVCRT')

"""
# ===================================================================== #
# ===================================================================== #

__all__ = ['scanf', 'sscanf', 'fscanf']

from __future__ import print_function
import re
import sys
from ctypes import (     # analysis:ignore
    util, CDLL, create_string_buffer, create_unicode_buffer, byref,
    c_byte, c_ubyte,
    c_short, c_ushort,
    c_int, c_uint,
    c_long, c_ulong,
    c_longlong, c_ulonglong,
    c_size_t,
    c_float, c_double, c_longdouble,
    c_char_p, c_wchar_p, c_void_p)

if sys.version_info < (3, 0):
    range = xrange     # analysis:ignore
else:
    unicode = str

LIBC = util.find_library('c')
if LIBC is None:
    LIBC = util.find_library('MSVCRT')  # windows version
# end if
LIBC = CDLL(LIBC)

C_SCANF_TYPES = {
    'i'   : c_int,
    'hhi' : c_byte,
    'hi'  : c_short,
    'li'  : c_long,
    'lli' : c_longlong,
    'ji'  : c_longlong,
    'zi'  : c_size_t,
    'ti'  : c_longlong,

    'd'   : c_int,
    'hhd' : c_byte,
    'hd'  : c_short,
    'ld'  : c_long,
    'lld' : c_longlong,
    'jd'  : c_longlong,
    'zd'  : c_size_t,
    'td'  : c_longlong,

    'u'   : c_uint,
    'hhu' : c_ubyte,
    'hu'  : c_ushort,
    'lu'  : c_ulong,
    'llu' : c_ulonglong,
    'ju'  : c_ulonglong,
    'zu'  : c_size_t,
    'tu'  : c_longlong,

    'o'   : c_uint,
    'hho' : c_ubyte,
    'ho'  : c_ushort,
    'lo'  : c_ulong,
    'llo' : c_ulonglong,
    'jo'  : c_ulonglong,
    'zo'  : c_size_t,
    'to'  : c_longlong,

    'x'   : c_uint,
    'hhx' : c_ubyte,
    'hx'  : c_ushort,
    'lx'  : c_ulong,
    'llx' : c_ulonglong,
    'jx'  : c_ulonglong,
    'zx'  : c_size_t,
    'tx'  : c_longlong,

    'f'   : c_float,
    'lf'  : c_double,
    'Lf'  : c_longdouble,
    'e'   : c_float,
    'le'  : c_double,
    'Le'  : c_longdouble,
    'g'   : c_float,
    'lg'  : c_double,
    'Lg'  : c_longdouble,
    'a'   : c_float,
    'la'  : c_double,
    'La'  : c_longdouble,

    'c'   : lambda l: lambda : create_string_buffer(l),  # c_char_p,
    'lc'  : lambda l: lambda : create_unicode_buffer(l),  # c_wchar_p,
    's'   : lambda l: lambda : create_string_buffer(l),  # c_char_p,
    'ls'  : lambda l: lambda : create_unicode_buffer(l),  # c_wchar_p,

    ']'   : lambda l: lambda : create_string_buffer(l),  # c_char_p,
    #'l[]' : c_wchar_p, handled in _get_c_object

    'p'   : c_void_p,

    'n'   : c_int,
    'hhn' : c_byte,
    'hn'  : c_short,
    'ln'  : c_long,
    'lln' : c_longlong,
    'jn'  : c_longlong,
    'zn'  : c_size_t,
    'tn'  : c_longlong,
}

SPECIFIER = re.compile('%([^ \t\n\r\f\v%%*]+)')


class IllegalSpecifier(Exception): pass


def _get_c_object(part, length):
    ctor = None

    # search most appropriate type constructor
    for pos in range(len(part)):
        try:
            ctor = C_SCANF_TYPES[part[-(pos+1):]]
        except KeyError:
            break
    if not ctor:
        raise IllegalSpecifier('cannot handle specifier "%%%s"' % part)

    # special handling of string types
    if part[-1:] in ('c', 's', ']'):
        # create unicode type for l[]
        if part[-1:] == ']' and part.find('l[') != -1:
            ctor = lambda l: lambda : create_unicode_buffer(l)
        # string buffers with length of input string
        ctor = ctor(length)
    return ctor()

# ===================================================================== #
# ===================================================================== #

def sscanf(s, fmt, outform=tuple):
    """
    clib sscanf for Python.
    For unicode strings use the l-versions of the string specifiers
    (%ls instead of %s).

    Returns a list with filled specifiers in order.
    """
    length = len(s)
    args = [_get_c_object(part, length) for part in SPECIFIER.findall(fmt)]
    sscanf_func = LIBC.sscanf    # scanf except: reads formatted input null-terminated string s
    buffer_ctor = create_string_buffer
    if isinstance(s, unicode):
        sscanf_func = LIBC.swscanf      # wscanf except: reads formatted input null-terminated string ws
        buffer_ctor = create_unicode_buffer
    filled = sscanf_func(buffer_ctor(s), buffer_ctor(fmt), *map(byref, args))
    return outform([args[i].value for i in range(filled)])

def scanf(fmt, outform=tuple):
    """
    clib scanf for Python.
    For unicode strings use the l-versions of the string specifiers
    (%ls instead of %s).

    Returns a list with filled specifiers in order.
    """
    return outform([sscanf(s, fmt, outform=list) for s in sys.stdin])

def fscanf(fid, fmt, size=None, outform=tuple):
    """
    clib fscanf for Python.
    """
    if size is None:
        return [sscanf(line, fmt) for line in fid]
    # end if
    return outform([sscanf(s, fmt, outform=list) for s in fid])

    # nrow = 1
    # ncol = 1
    # if len(_np.atleast_1d(size))>1:
#        nrow = size[0]
#        ncol = size[1]
    # else:
#        ncol = size
    # end if
#    return outform([[sscanf(s, fmt, outform=list) for jj in range(ncol)] for ii in range(nrow)])

#    s = sys.stdin
#    return sscanf(s, fmt)
#    if hasattr(s, "readline"): s = s.readline()
    # s = sys.stdin
#    length = len(s)    # this is problematic ... how big do I need to make str objects?
#    args = [_get_c_object(part, length) for part in SPECIFIER.findall(fmt)]
#    scanf_func = LIBC.scanf   # reads formatted input from stdin
#    buffer_ctor = create_string_buffer
#    if isinstance(s, unicode):
#        scanf_func = LIBC.wscanf
#        buffer_ctor = create_unicode_buffer
#    filled = scanf_func(buffer_ctor(fmt), *map(byref, args))
#    return tuple([args[i].value for i in range(filled)])

def ftell(fid):
    pass

def fgetl(fid):
    pass

def fgets(fid):
    pass

def fseek(fid):
    pass

# ===================================================================== #
# ===================================================================== #

if __name__ == '__main__':
    print(sscanf('abc defg %% bbb', '%s %s %%%% %s'))
    print(sscanf(u'abc defg %% bbb', u'%ls %ls %%%% %s'))

    print(sscanf(u'äüöß', u'%ls'))

    print(sscanf('ttttt abc - 123 -123.12345e-12 1b', '%5c %s - %d %f %x'))
    print(sscanf(u'ttttt abc - 123 -123.12345e-12 1b', u'%5lc %ls - %d %f %x'))

    print(sscanf('tttttabc', '%*5c%s'))
    print(sscanf(u'tttttabc', u'%*5lc%s'))

    print(sscanf(u'ääääääääääää 1', u'%3l[ä]%*l[ä] %d'))
# end if

# ===================================================================== #
# ===================================================================== #

