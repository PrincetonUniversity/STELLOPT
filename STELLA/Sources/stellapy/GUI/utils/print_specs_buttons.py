
def print_specs_buttons(style):
    background=style.lookup('TLabel', 'background')
    foreground=style.lookup('TLabel', 'foreground')
    activebackground=style.lookup('TLabel', 'background', ['active'])
    activeforeground= style.lookup('TLabel', 'foreground', ['active'])
    disabledforeground=style.lookup('TLabel', 'foreground', ['disabled'])
    print("Label",background,foreground,activebackground,activeforeground,disabledforeground)
    background=style.lookup('TButton', 'background')
    foreground=style.lookup('TButton', 'foreground')
    activebackground=style.lookup('TButton', 'background', ['active'])
    activeforeground= style.lookup('TButton', 'foreground', ['active'])
    disabledforeground=style.lookup('TButton', 'foreground', ['disabled'])
    borderwidth=style.lookup('TButton', 'borderwidth')
    # #232829 #ffffff #474d4e #ffffff #b3b3b3
    print("Button",background,foreground,activebackground,activeforeground,disabledforeground,borderwidth)


