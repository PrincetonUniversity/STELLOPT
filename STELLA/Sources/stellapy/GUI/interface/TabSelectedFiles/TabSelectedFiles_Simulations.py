####################################################################
#              CLASS FOR THE FRAME "SIMULATIONS"
####################################################################
''' SUBFRAME SIMULATIONS ON TAB 1

Manage the "Simulations" LabelFrame on the top left side of the tab "Simulations" 
which is attached to the root. Make it possible to select simulations and display 
the selection in the simulations_frame.


Absolute path of class Simulations
----------------------------------
root.tab_Simulations.class_simulations
     
'''

import os, pathlib
import tkinter as tk
from bisect import bisect
from tkinter import ttk
from stellapy.utils import get_filesInFolder
from stellapy.config import CONFIG, display_information #@unresolvedimports
from stellapy.GUI.utils import get_initialDirectory
from stellapy.GUI.interface import update_GUI
from stellapy.GUI.graph_tools.tkfilebrowser.functions import askopendirnames, askopenfilenames 

####################################################################
# INITIALIZE THE FRAME WITH INFORMATION ON THE SELECTED SIMULATIONS
####################################################################

class Simulations():
    '''
    Calling this class initiates the "Simulations" LabelFrame on the top left side 
    of the tab "Simulations" attached to the root.
    
       
    Attributes
    ----------
    input_files: dict[folder][input_file] = tk.IntVar()
        Stores the absolute path of the simulations, grouped by their folder, as well as wether or not to plot them.
           
    selected_files: list of str
        Saves which files were selected through select_files() and select_folders()
            

    Parent widgets
    --------------
    root:               Tk()
    tabheader:          ttk.Notebook(root)
    tab1:               ttk.Frame(tabheader)          
    frame_left:         ttk.Frame(tab1)
    frame_simulation:   ttk.LabelFrame(frame_left)    --> frame
    
      
    Widgets
    -------
    rightclick_menu: tk.Menu(root)
        Linked to the simulation labels, opens a menu which allows to open the simulation files
        
    button_selectFiles: ttk.Button(frame_simulation, select_files())
        Button to select simulations
        
    button_selectFolders: ttk.Button(frame_simulation, select_folders()) 
        Button to select folders which will select all simulations inside
        
    button_clearSimulations: ttk.Button(frame_simulation, clear_simulations())
        Button to clear all the selected simulations
    '''

#################################################################
#                          WIDGETS
#################################################################

    def __init__(self,tab1):
        '''
        Initialize the widgets in the "Simulations" LabelFrame of the tab "Simulations".
            - A scrollable canvas which will display the selected simulations
            - Thee buttons "Select simulations", "Select folders" and "Clear simulations"
            - A popup menu that shows up when right-clicking an input-file, this will
                give the possibility to open the corresponding data files.
                
        Initialize the class variables that will be used by the class methods.
            - input_files to save which simulations were selected
            - selected_files to save the newly selected simulations
        '''


        #====================================================================
        # Safe info from the root so that it can be passed on or used locally
        #====================================================================
        
        self.tab1 = tab1            
        self.root = tab1.root 
        self.frame = tab1.frame_simulation
                    
        #===========
        # VARIABLES
        #===========

        # Some local variables
        self.tree_folders = []
        self.selected_files = []                      
        self.initialdir = get_initialDirectory()   
        
        #============================
        # WIDGETS CREATION FOR FRAME
        #============================

        # The list of simulations will be displayed inside a treeview widget
        self.tree = ttk.Treeview(self.frame, show="tree",columns=("#0","#1","#2"))

        # "Select simulations" or "Add more simulations" button which updates the scrollable widget
        self.button_selectFiles = ttk.Button(self.frame, text="Select Simulations", width=20)
        self.button_selectFiles.config(command=lambda: self.select_files())

        # "Select Folder" or "Add more folders" button which updates the scrollable widget
        self.button_selectFolders = ttk.Button(self.frame, text="Select Folders", width=20)
        self.button_selectFolders.config(command=lambda: self.select_folders())

        # "Clear self" button which deletes the selected self and updates the scrollable widget
        self.button_clearSimulations = ttk.Button(self.frame, text="Clear simulations", width=20)
        self.button_clearSimulations.config(command=lambda: self.clear_simulations())
    
        
        #=================
        # CONFIGURE FRAME
        #=================
        tk.Grid.rowconfigure(   self.frame, 0, weight=1) # Scrollable canvas
        tk.Grid.rowconfigure(   self.frame, 1, weight=0) # 3 buttons to edit the simulation selection
 
        tk.Grid.columnconfigure(self.frame, 0, weight=1) # Column for button_selectFiles
        tk.Grid.columnconfigure(self.frame, 1, weight=1) # Column for button_selectFolders
        tk.Grid.columnconfigure(self.frame, 2, weight=1) # Column for button_clearSimulations

        #======================
        # WIDGETS ARRANGEMENT
        #=====================
        self.tree.grid(                      in_=self.frame, row=0, column=0, padx=8, pady=4, sticky='nesw', columnspan=3)
        self.button_selectFiles.grid(        in_=self.frame, row=1, column=0, padx=8, pady=4, sticky='nw', ipady=7)
        self.button_selectFolders.grid(      in_=self.frame, row=1, column=1, padx=8, pady=4, sticky='nw', ipady=7) 
        self.button_clearSimulations.grid(   in_=self.frame, row=1, column=2, padx=8, pady=4, sticky='nw', ipady=7) 
        
        #============
        # POPUP MENU
        #============
        self.rightclick_menu = tk.Menu(self.root, tearoff=0)
        self.rightclick_menu.add_command(label='Open ".in" file',            command = lambda: self.open_file(".in"))
        self.rightclick_menu.add_command(label='Open ".geometry" file',      command = lambda: self.open_file(".geometry"))
        self.rightclick_menu.add_command(label='Open ".fluxes" file',        command = lambda: self.open_file(".fluxes"))
        self.rightclick_menu.add_command(label='Open ".omega" file',         command = lambda: self.open_file(".omega"))
        self.rightclick_menu.add_command(label='Open ".final_fields" file',  command = lambda: self.open_file(".final_fields"))

        #============
        # TREE VIEW
        #============
            
        # Create a treeview widget and bind keypresses to the widget
        self.tree.bind("<Delete>", self.deleterow)
        self.tree.bind("<BackSpace>", self.deleterow)
        self.tree.bind("<Button-3>", self.popup_menu)
        
        # Add columns  
        self.tree.heading("#0", text="  Folders",anchor=tk.W)
        self.tree.heading("#1", text="Experiment", anchor=tk.W) 
        self.tree.heading("#2", text="Simulation", anchor=tk.W) 
        self.tree.heading("#3", text="Full path", anchor=tk.W) 
        self.tree["displaycolumns"] = ("#0", "#1", "#2")

        # When the tab is visible, make sure the columns have a good size
        def resize_columns(*args): 
            self.tree.column("#0", width=int(self.tree.winfo_width()*4/10), stretch=tk.YES)
            self.tree.column("#1", width=int(self.tree.winfo_width()*3/10), stretch=tk.YES)
            self.tree.column("#2", width=int(self.tree.winfo_width()*3/10), stretch=tk.YES) 
        self.frame.bind("<Visibility>", resize_columns)
    
        # Prevent folding of header comment on next line
        if True: return
                                   
#################################################################
#                          METHODS
#################################################################

#====================
# For the tree view
#====================

    def deleterow(self, event):
        '''Delete selected row'''
        
        # Get the information of the selected row
        if len(self.tree.selection()) != 0: 
            rowiid_parent = self.tree.selection()[0]
            rowiids = self.tree.get_children(item=rowiid_parent)
            
            # Delete all the child notes if the selected row was a folder
            for rowiid in rowiids:
                
                # Get the values of the row
                values = list(self.tree.item(rowiid, 'values'))
                
                # Delete the selected row
                try: self.tree.delete(rowiid)
                except: pass
                
                # Delete the corresponding input file
                try: self.root.input_files.remove(pathlib.Path(values[1]))
                except: pass
                
            # Now delete the parent node
            values = list(self.tree.item(rowiid_parent, 'values'))
            
            # Delete the input file or the folder
            try: # Delete the corresponding input file
                self.root.input_files.remove(pathlib.Path(values[1]))
            except: # Detele the folder
                index = None
                for tree_folder in self.tree_folders: 
                    self.tree.selection_set(tree_folder)  
                    folder_iid = self.tree.selection()[0]
                    if folder_iid == rowiid_parent:
                        index = self.tree_folders.index(tree_folder)
                self.tree_folders.pop(index)

            # Delete the selected row
            try: self.tree.delete(rowiid_parent)
            except: pass
                        
            # Update the GUI
            update_GUI(self.root)
        if True: return 
        
#====================
# For the popup menu
#====================

    def popup_menu(self, event):
        ''' Rightclick event linked to each self._widgets[folder]['label'][i]
        
        When rightclicking on a simulation label, it grabs its input_file and 
        folder and open up the popup menu which allows to open files.
        
        
        Attributes
        ----------
        rightclick_menu_input_file, rightclick_menu_folder: str
            Save the clicked [input_files, folder] so that the pop-up menu knows which file to open
        '''

        # Get the information of the selected row
        iid = self.tree.identify_row(event.y)
        
        # If the right click happened on a row: open the right click menu
        if iid:
            self.tree.selection_set(iid)     
            if len(self.tree.selection()) != 0:
                rowid = self.tree.selection()[0]
            values = list(self.tree.item(rowid, 'values'))
            
            # Save the selected input file
            self.rightclick_menu_input_file = pathlib.Path(values[1])
            
            # Pop-up the menu at the mouse location
            try:        self.rightclick_menu.tk_popup(event.x_root, event.y_root)
            finally:    self.rightclick_menu.grab_release()
        return 

    #-----------------------------
    def open_file(self, extension):
        ''' Command executed when choosing an option of self.rightclick_menu
        
        After rightclicking a simulation, the input file and folder were saved
        and a pop-up menu is opened, when choosing an option in the pop-up menu, 
        open the file with the corresponding extension, linked to the choosen input_file.
        '''
        
        # Get which simulation was right-clicked
        input_file = self.rightclick_menu_input_file
        
        # Construct the path to the file that needs to be opened
        file_name = input_file.with_suffix(extension)
        
        # If the file exists, open it with the texteditor written in the CONFIG file
        if os.path.isfile(file_name):
            if os.system(CONFIG['GENERAL']['texteditor'] + ' ' + str(file_name)) != 0: # Return code is non-zero if it failed to execute
                display_information(self.root,'Error','The text editor '+CONFIG['GENERAL']['texteditor'] + ' is not installed, chose another one.')
        if not os.path.isfile(file_name):
            display_information(self,'Warning','There is no "' + extension + '" file corresponding to ' + input_file + '.')

#====================
# Select simulations
#====================

    def select_files(self):
        ''' 
        Select simulations by choosing multiple files with the extension ".in".
        If simulations where selected, update the GUI through update_GUI().
        
        
        Attributes
        ----------
        selected_files: list of str
            Save the selected input files.
        '''

        # Choose files and start selection in the standard run folder.
        title           = "Select input file"
        filetypes       = (("in files","*.in"),("all files","*.*")) 
        selected_files  = askopenfilenames(initialdir=self.initialdir["Files"],title=title,filetypes=filetypes)

        # Update the GUI
        if selected_files is not None:
            if len(selected_files) != 0:
                selected_files  = [ pathlib.Path(i) for i in selected_files ]
                self.root.input_files = list(set(self.root.input_files + selected_files))
                update_GUI(self.root)
        return 
    
    #------------------------
    def select_folders(self):
        ''' 
        Select simulations by selecting one folder at a time and reading all input files inside of them.
        If simulations where selected, update the GUI through update_GUI().
        
        
        Attributes
        ----------
        selected_files: list of str
            Save the selected input files.
        '''

        # Choose files and start selection in the standard run folder.
        title             = "Select folder"
        selected_folders  = []
        selected_folders  = askopendirnames(title=title, initialdir=self.initialdir["Folders"])

        # Extract the files in the folder
        if selected_folders is not None:
            selected_folders  = [ pathlib.Path(f) for f in selected_folders ]
            selected_files = get_filesInFolder(selected_folders, end=".in")
    
            # Update the GUI
            if selected_files is not None:
                if len(selected_files) != 0:
                    self.root.input_files = list(set(self.root.input_files + selected_files))
                    update_GUI(self.root)   
        return 
    
    #--------------------
    def select_directories(self, title: str, selected_directory_list: list, initialdir):
        ''' Will reopen the directory window to add more folders. '''
        directory_path_string = tk.filedialog.askdirectory(initialdir=initialdir, title=title)
        if len(directory_path_string) > 0:
            selected_directory_list.append(directory_path_string)
            self.select_directories('Select the next Directory or Cancel to end', 
                                   selected_directory_list,
                                   os.path.dirname(directory_path_string))
            return selected_directory_list
        return 
    
    #---------------------------
    def clear_simulations(self):
        '''
        Remove all the selected simulations.
        Remove corresponding folder titles and input_files labels, crosses and checkmarks in the simulations labelFrame.
        Next update the GUI through update_GUI().
        
        
        Attributes
        ----------
        input_files: dict[folder][input_file] = tk.IntVar()
            Remove the input_files, fodlers and variables.
            
        widgets: dict['folder']['title', 'label', 'button', 'check'] = list of widgets; dict['folder']['var_lbl'] = list of tk.StringVar
            Destroy the widgets of each input_file, remove the widgets from the list and destroy the key in the dictionary.
            
        frame_folders: dict['folder'] = ttk.Frame(master=self.scrollableCanvas.scrollable_frame)
            Destroy the frames and remove them from the dictionary.
            
        '''
        
        # Remove all the items from the tree view
        self.tree.delete(*self.tree.get_children())
        self.tree_folders = []
        
        # Remove the labels from the research frame
        self.root.tab_Simulations.class_research.clear_simulations()
        
        # Remove the input files
        self.root.input_files = []
        
        # Update the GUI
        update_GUI(self.root)
        if True: return

#====================
# Update simulations
#====================

    def update_treeView(self):
        '''
        Update the scrollable canvas which displays the selected input files, 
        the selection is sorted by the parents folders. First make a frame
        and title label for each folder, next add the label, cross button and 
        checkmark widget for each input file or overwrite the existing labels
        with the correct data, overwriting reduces flickering of the GUI.
        
        Attributes
        ----------
        input_files: dict[folder][input_file] = tk.IntVar()
            Go thorough the simulations to display them on the GUI
            
        widgets : dict['folder']['title', 'label', 'button', 'check'] = list of widgets; dict['folder']['var_lbl'] = list of tk.StringVar
            Stores the newly created widgets for the simulations 
            
        frame_folders: dict['folder'] = ttk.Frame(master=self.scrollableCanvas.scrollable_frame)
            Stores the newly created widgets for the folders
        '''

        # Only update when there are input files
        if len(self.root.input_files) != 0:
    
            # Sort simulations by folder
            unique_folders = self.root.Research.unique_folders
    
            # Iterate over the folders
            for folder in unique_folders:
                
                # Look at the files in the current <folder>
                files = []; input_files = []; folders = []; e_id = []; s_id = []
                for i in range(len(self.root.Research.files)):
                    if self.root.Research.folders[i]==folder:
                        files.append(self.root.Research.files[i])
                        input_files.append(self.root.Research.input_files[i])
                        folders.append(self.root.Research.folders[i])
                        try: 
                            e_id.append(self.root.Research.e_id[i])
                            s_id.append(self.root.Research.s_id[i])
                        except: 
                            e_id.append("New GUI feature")
                            s_id.append("New GUI feature")
    
                # Check whether the folder is already in the treeview
                # If it is not, add the folder to the treeview
                if not self.tree.exists(folder):
                    contents = [self.tree.item(child)["text"] for child in self.tree.get_children("")]
                    self.tree_folders.append(self.tree.insert(parent="", index=bisect(contents, folder), iid=folder, text=folder, values=("","")))
        
                # Add the input files to the folder
                for tree_folder in self.tree_folders:
                    self.tree.selection_set(tree_folder)  
                    folder_iid = self.tree.selection()[0]
                    if folder_iid==folder:
                        for input_file in input_files:
                            # If the input file is not in the treeview, add it
                            if not self.tree.exists(str(input_file)):
                                i = input_files.index(input_file)
                                contents = [self.tree.item(child)["text"] for child in self.tree.get_children(folder)]
                                self.tree.insert(tree_folder, bisect(contents, files[i]), iid=str(input_file), text=files[i], values=(s_id[i], e_id[i], str(input_file)))           
                            # If the input file is already in the treeview, make sure it has the correct columns
                            if self.tree.exists(str(input_file)):
                                i = input_files.index(input_file)
                                self.tree.item(str(input_file), text=files[i], values=(s_id[i], e_id[i], str(input_file)))
             
            # Remove focus from items
            for item in self.tree.selection():
                self.tree.selection_remove(item)



