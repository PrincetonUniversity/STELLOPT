


def print_decoratorOnCommandPrompt(name, style="long"):
    ''' When a function is called print ---- \\n <function_name> \\n ---- ''' 

    # If a function is called within a function, indent the prints in the terminal
    indent_string = "     "

    # Header to be printed in the terminal when a function is called
    if name == "begin":

        # indent the terminal prints
        indent(1)

        if style == "long":
            width = 60
            print(indent_string*indent.counter, "".center(width,"-"))

        if style == "short":
            pass

    # Ending line to indicate the end of a funtion
    elif name == "end":

        if style == "long":
            width = 60
            print(indent_string*indent.counter, "".center(width,"-"))
            print()

        if style == "short":
            pass 
        
        # unindent the terminal prints
        indent(-1)

    # Header with module.function to be printed in the terminal when a function is called
    else:

        # indent the terminal prints
        indent(1)

        if style == "long":
            width = 60
            print()
            print(indent_string*indent.counter, "".center(width,"="))
            print(indent_string*indent.counter, name.center(width))
            print(indent_string*indent.counter, "".center(width,"="),"\n")

        if style == "short":
            width = len(name);
            print("    "*indent.counter, "".center(width,"-"), "\n", name.center(width), "\n", "".center(width,"-"), sep="")

    return indent_string*indent.counter


def indent(amount): 
    try:    indent.counter += amount;
    except: indent.counter = 0;

 

