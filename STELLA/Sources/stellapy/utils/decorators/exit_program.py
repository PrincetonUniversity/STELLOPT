
import sys

def exit_program(reason, function, code_line):
    ''' Exit the program and print the reason to the command prompt.
    
    Example
    -------
    exit_program(reason, function_name, sys._getframe().f_lineno) 
     '''

    from stellapy.utils.decorators.verbose_wrapper import indent

    width = 60
    print()
    print(indent, 'EXIT PROGRAM   -->  ', reason); 
    print(indent, '               -->  ', function.__module__);
    print(indent, '               -->   line', code_line);
    print(indent, "".center(width,"-"))
    sys.exit(2)
