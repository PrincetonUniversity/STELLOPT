
#################################################################
#                   CLASS FOR THE FIRST TAB
#################################################################
''' 
'''

# Load modules
import os
import tkinter as tk
import numpy as np
from tkinter import ttk
from tkinter import messagebox

# Personal modules

#################################################################
#                   CLASS FOR THE FIRST TAB
#################################################################
class TabVideos:
 
#################################################################
#                          WIDGETS
#################################################################

    # Initate the tab for "Select simulations"
    def __init__(self, tabMaster): 


        #======================================================
        # Safe info from the tab and root that can be passed on
        #======================================================

        self.window = tabMaster                      # Direct accesses to the window of tab "Select Simulations" 
        self.root = tabMaster.root                   # Needed to center windows based on the root screen

        #====================
        # Create the subtabs
        #====================

        # Make multiple tabs
        tab_control = ttk.Notebook(self.window, style='subheader.TNotebook')
        tab1 = ttk.Frame(tab_control) 
        tab2 = ttk.Frame(tab_control) 
        tab3 = ttk.Frame(tab_control)  
        tab_control.add(tab1, text='Time evolution of the Potential')  
        tab_control.add(tab2, text='Time evolution of omega/gamma') 
        tab_control.add(tab3, text='Time evolution the fluxes')  
        tab_control.pack(expand=1, fill='both')

        # Create the interfaces
        for tab in [tab1, tab2, tab3]:
            tab.root = self.root # Make the root accessible for the tabs
            
        # Add treeview to tab 1
        self.tree = tree = ttk.Treeview(tab1)
        
        # Bind the selection of a cell
        self.tree.bind("<<TreeviewSelect>>", self.selectItemOnClick)
        self.tree.bind("<Double-1>", self.onDoubleClick)
        self.tree.bind("<Delete>", self.deleterow)
        self.tree.bind("<BackSpace>", self.deleterow)
        
        # Add columns
        tree["columns"]=("one","two","three")
        tree.column("#0",    width=270, minwidth=270, stretch=tk.YES)
        tree.column("one",   width=150, minwidth=150, stretch=tk.YES)
        tree.column("two",   width=400, minwidth=200, stretch=tk.YES)
        tree.column("three", width=80,  minwidth=50,  stretch=tk.YES)
        
        # Add heading
        tree.heading("#0", text="Name",anchor=tk.W)
        tree.heading("#1", text="Date modified",anchor=tk.W)
        tree.heading("#2", text="Type",anchor=tk.W)
        tree.heading("#3", text="Size",anchor=tk.W)
        
        # Insert information: tree.insert(parent, index, iid, text, values)
        # If you want the parent widget as the master (root) node, we can set this to the empty string (”)
        folder1 = tree.insert(parent="", index=1, iid="folder1", text="Folder 1", values=("23-Jun-17 11:05","File folder",""))
        tree.insert("", 2, text="text_file.txt", values=("23-Jun-17 11:25","TXT file","1 KB"))
        # Level 2
        for i in range(50):
            tree.insert(folder1, "end", text="photo"+str(i)+".png", values=("23-Jun-17 11:28","PNG file",str(i)+" KB"))
#         tree.insert(folder1, "end", text="photo2.png", values=("23-Jun-17 11:29","PNG file","3.2 KB"))
#         tree.insert(folder1, "end", text="photo3.png", values=("23-Jun-17 11:30","PNG file","3.1 KB"))
            
        # Add the widget to the tab
        #tree.pack(side=tk.TOP,fill=tk.X)

    def deleterow(self, event):
        '''Delete selected row'''
        
        # Close the popup if there is one
        try: self.entryPopup.on_return()
        except: pass
        
        # Now delete the selected row
        if len(self.tree.selection()) != 0:
            row = self.tree.selection()[0]
        try: self.tree.delete(row)
        except: pass
                
    def selectItemOnClick(self, *args):
        ''' Get the value in column 0 of the selected row. '''
        selitems = self.tree.selection()
        if selitems:
            selitem = selitems[0]
            text = self.tree.item(selitem, "text")  

    def onDoubleClick(self, event):
        ''' Executed, when a row is double-clicked. Opens EntryPopup above the item's column, 
        so it is possible to change the text. '''
    
        # Close previous popup if there is one
        try: self.entryPopup.on_return()
        except: pass
    
        # Select row and column that was clicked on
        rowid = self.tree.identify_row(event.y)
        columnid = self.tree.identify_column(event.x)
    
        # Get column position info
        x,y,width,height = self.tree.bbox(rowid, columnid)
    
        # Place Entry popup properly         
        text = self.tree.item(rowid, 'text')
        values = self.tree.item(rowid, 'values')
        self.entryPopup = EntryPopup(self.root, self.tree, rowid, columnid, text, values, self.update_tree)
        self.entryPopup.place(x=x, y=y+height//2, anchor=tk.W)
        
    def update_tree(self):
        print("DO SOMETHING")

 
class EntryPopup(ttk.Entry):

    def __init__(self, root, tree, iid, columnid, text, values, update_tree, **kw):
        super().__init__(tree, style='tree.TEntry', **kw)
        
        # Save the information about the treeview
        self.tree = tree
        self.update_tree = update_tree
        self.iid = iid
        self.columnid = columnid
        self.text = text
        self.values = list(values)
        self['exportselection'] = False
        
        # Insert the text of the cell in the entry widget: either its text or one of the values
        if "0" in self.columnid: self.insert(0, text) 
        for i in range(1,len(values)+1):
            if str(i) in self.columnid: self.insert(0, values[i-1]) 

        # Bind the key presses to the entry widget
        self.focus_force()
        self.bind("<Return>", self.on_return)
        self.bind("<Control-a>", self.select_all)
        self.bind("<Escape>", lambda *ignore: self.destroy())

    def on_return(self, event=None):
        if "0" in self.columnid: self.text = self.get()
        for i in range(1,len(self.values)+1):
            if str(i) in self.columnid: 
                self.values[i-1] = self.get()
        self.tree.item(self.iid, text=self.text, values = self.values) 
        self.update_tree()
        self.destroy()

    def select_all(self, *ignore):
        ''' Set selection on the whole text '''
        self.selection_range(0, 'end')

        # returns 'break' to interrupt default key-bindings
        return 'break'




