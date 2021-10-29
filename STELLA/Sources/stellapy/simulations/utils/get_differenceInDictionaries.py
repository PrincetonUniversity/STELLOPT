
from stellapy.utils import initiate_nesteddict

def get_differenceInDictionaries(dict1, dict2):
    ''' Returns the difference in keys between two dictionaries. '''
    
    # Initiate the different dictionary
    dict_difference = initiate_nesteddict()
    
    # Get the difference of a two layered dictionary
    if isinstance(dict1[list(dict1.keys())[0]], dict):
        for knob, sub_dict1 in dict1.items():
            for key, value in sub_dict1.items():
                if dict2[knob][key] != value:
                    if knob not in dict_difference:
                        dict_difference[knob] = [key]
                    else:
                        dict_difference[knob].append(key)
    
    return dict_difference
