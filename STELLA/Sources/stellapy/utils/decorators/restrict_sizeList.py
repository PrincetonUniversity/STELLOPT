#===============================================================
# Restrict size of printed list
#===============================================================

import numpy as np

def restrict_sizeList(input_list, size):
    ''' Print a long list as [ a b c ... x y z ] '''

    # For 1D arrays or list, return [ a b c ... x y z ]
    if np.array(input_list).ndim == 1:

        # Make sure it is a list
        input_list = list(input_list)

        # If the length is bigger than size*2, write it as [ a b c ... x y z ]
        if len(input_list) > size*2:
            
            string = "["
            for i in np.arange(size):
                string += "{0:.2f}".format(input_list[i]) + ", "
            string += "..., "
            for i in np.arange(len(input_list)-size, len(input_list)):
                string += "{0:.2f}".format(input_list[i]) + ", "
            return string[:-2] + "]"

        if len(input_list) == 1:
            return str(input_list[0])

        else:
            return str(input_list)


    # For 3D array, return [[[a, b, c], ... , [d, e, f]]]
    if np.array(input_list).ndim == 3:
        
        # Make sure it is an array
        input_list = np.array(input_list)
        string = np.array2string(input_list, formatter={'float_kind':lambda x: "%.2E" % x})
        return str(string[0:43]).replace("\n","") + ', ..., ' + str(string[-44:]).replace("\n","")






