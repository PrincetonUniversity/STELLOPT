
def print_inTwoColumns(name1, value1, unit1=None, name2=None, value2=None, unit2=None):
    ''' Print two columns of <name> <value> <unit> '''

    from stellapy.utils.decorators.verbose_wrapper import indent
    
    # Space between {name value} and between the two columns; format for the floats
    spacing1 = '{0:13}  {1:8}  {2:10}'
    spacing2 = '{0:10}'.format(" ")
    spacing3 = '{0:13}  {1:8}  {2:6}'
    float_format = "{:8.4f}"
    sci_format = "{:8.2E}"

    # If the values are floats, give them the correct format
    if isinstance(value1, float):
        if (value1 < 10**3):
            value1 = float_format.format(value1)
        elif (value1 >= 10**3):
            value1 = sci_format.format(value1)
    if isinstance(value2, float):
        if value2 < 10**3:
            value2 = float_format.format(value2)
        elif value2 >= 10**3:
            value2 = sci_format.format(value2)

    # If a title is given, value="" and name can be split over name and value
    if not unit1 and not name2 and not value2 and not unit2:
        print(indent, '{0:35}'.format(name1) + spacing2 + '{0:35}'.format(value1))
    else:
        print(indent, "    ", spacing1.format(name1, value1, unit1) + spacing2 + spacing3.format(name2, value2, unit2))



def print_inTwoColumns2(name1, value1, name2=None, value2=None):
    ''' Print two columns of <name> <dim> ; <name> <vector>'''

    from stellapy.utils.decorators.verbose_wrapper import indent
    
    # Space between {name value} and between the two columns; format for the floats
    spacing1 = '{0:20}  {1:12}'
    spacing2 = '{0:2}'.format(" ")
    spacing3 = '{0:2}  {1:30}'
    float_format = "{:8.4f}"
    sci_format = "{:8.2E}"

    # If the values are floats, give them the correct format
    if isinstance(value1, float):
        if (value1 < 10**3):
            value1 = float_format.format(value1)
        elif (value1 >= 10**3):
            value1 = sci_format.format(value1)
    elif value1:
        value1 = str(value1)

    if isinstance(value2, float):
        if value2 < 10**3:
            value2 = float_format.format(value2)
        elif value2 >= 10**3:
            value2 = sci_format.format(value2)    
    elif value2:
        value2 = str(value2)

    # If a title is given, value="" and name can be split over name and value
    if not value2:
        print(indent, '{0:20}'.format(name1) + spacing2 + '{0:15}'.format(value1) + spacing2 + '{0:42}'.format(name2))
    else:
        print(indent, "  ", spacing1.format(name1, value1) + spacing2 + spacing3.format(name2, value2))





def print_inTwoColumns3(name1, value1, name2=None, value2=None):
    ''' Print two columns of <name> <value> ; <name> <value>'''

    from stellapy.utils.decorators.verbose_wrapper import indent
    
    # Space between {name value} and between the two columns; format for the floats
    spacing1 = '{0:15}  {1:<8}'
    spacing2 = '{0:10}'.format(" ")
    spacing3 = '{0:15}  {1:<8}'
    float_format = "{:<8.4f}"
    sci_format = "{:<8.2E}"

    # If the values are floats, give them the correct format
    if isinstance(value1, float):
        if (value1 < 10**3):
            value1 = float_format.format(value1)
        elif (value1 >= 10**3):
            value1 = sci_format.format(value1)
    elif value1:
        value1 = str(value1)

    if isinstance(value2, float):
        if value2 < 10**3:
            value2 = float_format.format(value2)
        elif value2 >= 10**3:
            value2 = sci_format.format(value2)    
    elif value2:
        value2 = str(value2)

    # If a title is given, value="" and name can be split over name and value
    if not value2:
        print(indent, '{0:20}'.format(name1) + spacing2 + '{0:15}'.format(value1) + spacing2 + '{0:42}'.format(name2))
    else:
        print(indent, "  ", spacing1.format(name1, value1) + spacing2 + spacing3.format(name2, value2))
