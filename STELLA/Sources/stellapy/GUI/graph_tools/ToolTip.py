import tkinter as tk

class ToolTip(object):
    """
    create a tooltip for a given widget
    """
    def __init__(self, root, widget, text='widget info'):
        self.root = root
        self.waittime = 500     #miliseconds
        self.wraplength = 180   #pixels
        self.widget = widget
        self.text = text
        self.widget.bind("<Enter>", self.enter)
        self.widget.bind("<Leave>", self.leave)
        self.widget.bind("<ButtonPress>", self.leave)
        self.id = None
        self.tw = None

    def enter(self, event=None):
        self.schedule()

    def leave(self, event=None):
        self.unschedule()
        self.hidetip()

    def schedule(self):
        self.unschedule()
        self.id = self.widget.after(self.waittime, self.showtip)

    def unschedule(self):
        id_ = self.id
        self.id = None
        if id_:
            self.widget.after_cancel(id_)

    def showtip(self, event=None):
        root = self.root
        x, y, cx, cy = self.widget.bbox("insert") #@unusedvariable
        x += self.widget.winfo_rootx() + 25
        y += self.widget.winfo_rooty() + 20
        # creates a toplevel window
        self.tw = tk.Toplevel(self.widget)
        # Leaves only the label and removes the app window
        self.tw.wm_overrideredirect(True)
        self.tw.wm_geometry("+%d+%d" % (x, y))
        label = tk.Label(self.tw, text=self.text, justify=tk.LEFT, background=root.color['bbg'], foreground=root.color['fg'],\
                        activebackground=root.color['bbg'],activeforeground=root.color['bbg'],disabledforeground=root.color['bbg'],
                        relief=tk.FLAT, borderwidth=2, highlightcolor="red", font=("TkDefaultFont", "10", "normal"))
        label.pack(ipadx=5, ipady=5)

    def hidetip(self):
        tw = self.tw
        self.tw= None
        if tw:
            tw.destroy()
    
 