

#################################################################
#                   PREFERENCE WINDOW
#################################################################

# Load modules
import tkinter as tk
from tkinter import ttk

# Load personal modules
from stellapy.config import write_configurationFile
from stellapy.GUI.interface.PreferencesWindow.TabAbout import TabAbout
from stellapy.GUI.interface.PreferencesWindow.TabAppearance import TabAppearance
from stellapy.GUI.interface.PreferencesWindow.TabGeneral import TabGeneral

#================
# MENU CREATION 
#================
class PreferenceWindow:

    def __init__(self, root):
        
        
        #===========
        # VARIABLES
        #===========
        
        # Attach the root so we can carry it into the functions, get the tab for its dimensions
        self.root = root
        tab1 = root.tab1
        tab1.update()
        
        # Get the width and height of the root window + title bars and of simply the root window
        self.height = root.winfo_height() # Height of the application minus the application title
        self.width =  root.winfo_width()  # Width of the application  
        
        window_height = tab1.winfo_height()                     # Height of the window minus the tab header
        window_width  =  tab1.winfo_width()                     # Width of the window
        outerFrame_width  = tab1.winfo_rootx() - tab1.winfo_x() # Top left x coordinate of the window excluding the outer-frame and including it
        outerFrame_height = tab1.winfo_rooty() - tab1.winfo_y() # Top y coordinate of the window excluding the outer-frame and including it
        header_width = self.width - window_width                # Pixels between root window and notebook
        header_height = self.height - window_height             # Height of the header of the notebook
        header_height = header_height -  header_width           # Height of the header of the notebook
        
        #===========
        # WIDGETS
        #===========
              
        # Pack the button in a frame so it's the same size as the notebook header
        self.frame_openPreferences = ttk.Frame(self.root, height=header_height, width=100, style="dot.TFrame")
        self.frame_openPreferences.pack_propagate(0) # The size of the frame controls the size of the button rather than visa versa
        self.frame_openPreferences.place(relx=1, rely=0, anchor="ne")
        
        # Create a button that will open the preferences window
        self.btn_openPreferences = ttk.Button(master=self.frame_openPreferences, text="  ...  ", style='dot.TButton')
        self.btn_openPreferences.config(command=lambda: self.open_preferencesWindow())
        self.btn_openPreferences.pack(expand=1, fill=tk.BOTH)
        

        
#==============================
# Open the preferences window
#==============================
    
    def open_preferencesWindow(self):
        
        # Create the preferences window 
        self.window_preferences = tk.Toplevel(self.root)
        self.window_preferences.title("Preferences") 
        
        # Center the new window in the screen
        winx = 500;  x = self.width/2  - winx/2
        winy = 500; y = self.height/2 - winy/2
        self.window_preferences.minsize(winx, winy)
        self.window_preferences.geometry("+%d+%d" % (x, y))
        
        # Create a tabbed view with the possible settings
        self.tab_header = ttk.Notebook(self.window_preferences, style='header.TNotebook')
        
        # Add frames to the tab_header which are the tab windows
        self.tab_general     = ttk.Frame(self.tab_header)  # Text editor
        self.tab_about       = ttk.Frame(self.tab_header)  # About me
        self.tab_appearance  = ttk.Frame(self.tab_header)  # Fonts and colors 
        
        # Add the tabs to the tab header
        self.tab_header.add(self.tab_general, text='General')
        self.tab_header.add(self.tab_appearance, text='Appearance')
        self.tab_header.add(self.tab_about, text='About stellapy')
        self.tab_header.pack(expand=1, fill='both')
        
        # Attach the root so the classes can acces them
        self.tab_general.root = self.root
        self.tab_about.root = self.root
        self.tab_appearance.root = self.root
        
        # Fill the tabs with widgets through classes
        self.tabGeneral     = TabGeneral(self.tab_general)
        self.TabAbout       = TabAbout(self.tab_about)
        self.tabAppearance  = TabAppearance(self.tab_appearance)
    
        # Closing event: apply the changes, then close the window
        def on_closing():
            
            # Apply the changes
            self.tabAppearance.apply_changesTheme()
            self.tabGeneral.apply_changesTextEdiror()
            self.tabGeneral.apply_changesPaths()
            
            # Destroy the tkinter window
            self.window_preferences.destroy()
        
        self.window_preferences.protocol("WM_DELETE_WINDOW", on_closing)
        





 


