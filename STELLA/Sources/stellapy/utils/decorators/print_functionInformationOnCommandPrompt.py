

import pathlib
from .print_decoratorOnCommandPrompt import print_decoratorOnCommandPrompt

def print_functionInformationOnCommandPrompt(function, default_arguments, arguments):
    ''' When a bash command is called, print the function name, its arguments and its changed arguments.''' 

    # Write that we are using a bash command
    print_decoratorOnCommandPrompt("BASH COMMAND") 
    
    # Print the function name and its default arguments
    print(" The following python function will be executed: ")
    print("     ", function)
    print("\n The function has the following default arguments: ")
    for name, value in default_arguments.items(): 
        if isinstance(value, str): value="'"+value+"'"
        if isinstance(value, list): 
            value = [str(i) for i in value]
            value="["+', '.join(value)+"]"
        if value==None:  value="None"
        if value==True:  value="True"
        if value==False: value="False"
        print('     ', '{0:20}'.format(name+":") + '{:<40}'.format(value))
    
    # If some arguments were defined through the bash command, show them
    if arguments is not None: 
        print("\n The arguments changed through the command prompt are:")
        for name, value in arguments.items():  
            if arguments[name] != default_arguments[name]:
                if isinstance(value, pathlib.PurePath): value=str(value)
                if isinstance(value, str): value="'"+value+"'"
                if isinstance(value, list): 
                    value = [str(i) for i in value]
                    value="["+', '.join(value)+"]"
                if value==None:  value="None"
                if value==True:  value="True"
                if value==False: value="False"
                print('     ', '{0:20}'.format(name+":") + '{:<40}'.format(value))
            
    # If the function is not called from --help then add a header at the end
    if arguments is not None:  print_decoratorOnCommandPrompt("end")
    

 

