# Load modules
import tkinter as tk 
from stellapy.GUI.utils.get_centerCoordinatesOfMonitor import get_centerCoordinatesOfMonitor

# Display a window with some information
def display_information(root,title=None,message=None):
    ''' Display a window with information.

    Parameters
    ----------
    self : class object
        self.root should be linked to the root window in order to center the window on the screen.
    message : str
        Message to be displayed.
    '''

    if not title:   title="Please insert a title"
    if not message: message="Please insert a message."
    information_window = tk.Toplevel(bg='#474d4e')
    information_window.title(title)
    root.eval(f'tk::PlaceWindow {str(information_window)} center')
    return

# Display a window with some information
def display_warning(root, message):
    ''' Display a window with a warning. '''
    if root:
        # Create a warning message in a new popup window
        information_window = tk.Toplevel(bg=root.color['bg'])
        information_window.title("Warning")
        label = tk.Label(information_window, text=message, bg=root.color['bg'], font=('arial', 11, 'bold'))
        label.pack(ipadx=50, ipady=10, fill='both', expand=True)
        # Make sure the new window is centered
        information_window.withdraw()
        information_window.update_idletasks()   
        x, y = get_centerCoordinatesOfMonitor(root, information_window) 
        information_window.geometry("+%d+%d" % (x, y))
        information_window.deiconify()
    return
 
 
