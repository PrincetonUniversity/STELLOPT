# -*- coding: utf-8 -*-
"""
Created on Tue May 17 23:23:29 2016

    This is a data object that allows seemless conversion between class objects
    and python dictionary formats

    To create a data structure (class object) from a python dictionary use the
    initialization method of Struct.

    To create a python dictionary from a data structure (class object) use the
    Struct method dict_from_class
        ex usage//
            data = {'x': x, 'y': y}
            DataObject = Struct(data)
            x = DataObject.x
            y = DataObject.y
            data =  DataObject.dict_from_class()
            x = data['x']
            y = data['y']

@author: Gavin Weir <gavin.weir@ipp.mpg.de>
"""
# ===================================================================== #
# ===================================================================== #

from __future__ import absolute_import, with_statement, absolute_import, \
                       division, print_function, unicode_literals

import numpy as _np
# import OrderedDict
__metaclass__ = type

# ===================================================================== #

#Convert a dictionary to an object structure and vice-versa
class Struct(object):
    def __init__(self, d={}):
        """Convert a dictionary to a class

        @param : entries Dictionary
        """
        if isinstance(d, (list, tuple)):
            setattr(self, "items",
                    [Struct(x) if isinstance(x, dict) else x for x in d])
        else:
            for a, b in d.items():
                if isinstance(b, dict):
                    setattr(self, a, Struct(b))
                elif isinstance(b, (list, tuple)):
                    setattr(self, a,
                            [Struct(x) if isinstance(x, dict) else x for x in b])
                else:
                    setattr(self, a, b)
#        if entries.__class__.__name__!='list':
#            self.__dict2vars__(**entries)
        # end if
    #end def __init__


    def __getitem__(self, key):
        return self.items[key]

    # ====================== #

    def __dict2vars__(self, **entries):
        """Inherit variables from a dictionary

        @param : entries Dictionary
        """
#        for ent in entries:
#            self.__dict__.update()
        self.__dict__.update(entries)
    #end def __dict2vars__

    # ====================== #


    def dict_from_class(self, excluded_keys = None):
        return self.__class2dict__(self, excluded_keys)
    #end def __dict_from_class__

    # ====================== #

    @staticmethod
    def __class2dict__(cls, _excluded_keys = None):
        """Convert the fields of this object/class to a dictionary

        @param :class:Struct
        @return :dictionary
        """
        if _excluded_keys is None:
#            _excluded_keys = Struct().__dict__.keys()
#            _excluded_keys = set(Struct().__dict__.keys())
#        else:
            _excluded_keys = ['None']
        # end if
#        tuples = []
#        for (key, value) in cls.__dict__.items():
#            if key not in _excluded_keys:
#                if type(value)!=type(Struct()):
#                    tuples.append((key, value))
#                else:
#                    tuples.append((key, value.dict_from_class()))
#                # end if
#            # end if
#        # end if
#        return dict(tuples)   # this is super inefficient

        def _check_excluded_keys_and_type(key, value):
            if key not in _excluded_keys and type(value)!=type(Struct()):
                return (key,value)
            elif key not in _excluded_keys:
                return (key, value.dict_from_class())
            # end if
        # end def
        return dict(_check_excluded_keys_and_type(key,value) for (key, value) in cls.__dict__.items())
#
#        return dict( (key,value) if  type(value)!=type(Struct())
#                else (key,value.dict_from_class()) for (key, value) in cls.__dict__.items())
##
#        return dict( (key, value)
#                    for (key, value) in cls.__dict__.items()
#                        if key not in _excluded_keys )
    #end def __dict_from_class__

    # ====================== #


    def __class2NPdict__(self, obj, _excluded_keys = None):
        """Convert the fields of this object/class to a dictionary

        @param :class:Struct
        @return :dictionary
        """
        if _excluded_keys is None:
#            _excluded_keys = set(Struct().__dict__.keys())
#        else:
            _excluded_keys = ['None']
        # end if

        return dict((key, self.__list2np__(value))
                        for (key, value) in obj.__dict__.items()
                            if key not in _excluded_keys)
    #end def __class2NPdict__

    # ====================== #


    @staticmethod
    def get_object(adict):
        """Convert a dictionary to a class

        @param :adict Dictionary
        @return :class:Struct
        """
        return Struct(adict)
    #end def get_object

    # ====================== #


    @staticmethod
    def __list2np__(value):
        if value.__class__.__name__=='list':
            dtype = type(value[0]).__name__
            if dtype.find('int') > -1:
                return _np.asarray(value, dtype=_np.int64, order='C')
            elif dtype.find('float') > -1:
                return _np.asarray(value, dtype=_np.float64, order='C')
            # end if
        # end if
        return value
    # end def __list2np__

    # ====================== #

#    def __repr__(self):
##        # python27
##        return '<%s>' % str('\n '.join('%s : %s' % (k, repr(v)) for (k, v) in self.__dict__.iteritems()))
#        # python3
#        return '<%s>' % str('\n '.join('%s : %s' % (k, repr(v)) for (k, v) in self.__dict__.items()))
#    #end def __repr__

    # ====================== #



    # ====================== #

#end class Struct

# ===================================================================== #
# ===================================================================== #
