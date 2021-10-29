# *** coding: utf-8 -*-
"""
tkfilebrowser - Alternative to filedialog for Tkinter
Copyright 2017 Juliette Monsel <j_4321@protonmail.com>

tkfilebrowser is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

tkfilebrowser is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


Tooltip and TooltipTreeWrapper classes to display the full path of a shortcut
when the mouse stays over long enough
"""


from .constants import tk, ttk
from sys import platform


class Tooltip(tk.Toplevel):
    """Tooltip to display when the mouse stays long enough on an item."""
    def __init__(self, parent, **kwargs):
        """
        Create Tooltip.

        Options:
            * parent: parent window
            * text: text (str) to display in the tooltip
            * compound: relative orientation of the graphic relative to the text
            * alpha: opacity of the tooltip (0 for transparent, 1 for opaque),
                     the text is affected too, so 0 would mean an invisible tooltip
        """
        tk.Toplevel.__init__(self, parent)
        self.transient(parent)
        if platform.startswith('linux'):
            self.attributes('-type', 'tooltip')
        self.attributes('-alpha', kwargs.get('alpha', 0.8))
        self.overrideredirect(True)
        style = kwargs.get('style', 'tooltip.tkfilebrowser.TLabel')

        bg = ttk.Style(self).lookup(style, 'background')
        self.configure(background=bg)

        self.label = ttk.Label(self, text=kwargs.get('text', ''),
                               style=style, compound=kwargs.get('compound', 'left'),
                               padding=kwargs.get('padding', 4))
        self.label.pack()

    def configure(self, **kwargs):
        if 'text' in kwargs:
            self.label.configure(text=kwargs.pop('text'))
        if 'image' in kwargs:
            self.label.configure(image=kwargs.pop('image'))
        if 'alpha' in kwargs:
            self.attributes('-alpha', kwargs.pop('alpha'))
        tk.Toplevel.configure(self, **kwargs)


class TooltipTreeWrapper:
    """Tooltip wrapper for a Treeview."""
    def __init__(self, tree, delay=1500, **kwargs):
        """
        Create a Tooltip wrapper for the Treeview tree.

        This wrapper enables the creation of tooltips for tree's items with all
        the bindings to make them appear/disappear.

        Options:
            * tree: wrapped Treeview
            * delay: hover delay before displaying the tooltip (ms)
            * all keyword arguments of a Tooltip
        """
        self.tree = tree
        self.delay = delay
        self._timer_id = ''
        self.tooltip_text = {}
        self.tooltip = Tooltip(tree, **kwargs)
        self.tooltip.withdraw()
        self.current_item = None

        self.tree.bind('<Motion>', self._on_motion)
        self.tree.bind('<Leave>', self._on_leave)

    def _on_leave(self, event):
        try:
            self.tree.after_cancel(self._timer_id)
        except ValueError:
            # nothing to cancel
            pass

    def add_tooltip(self, item, text):
        """Add a tooltip with given text to the item."""
        self.tooltip_text[item] = text

    def _on_motion(self, event):
        """Withdraw tooltip on mouse motion and cancel its appearance."""
        if self.tooltip.winfo_ismapped():
            x, y = self.tree.winfo_pointerxy()
            if self.tree.winfo_containing(x, y) != self.tooltip:
                if self.tree.identify_row(y - self.tree.winfo_rooty()):
                    self.tooltip.withdraw()
                    self.current_item = None
        else:
            try:
                self.tree.after_cancel(self._timer_id)
            except ValueError:
                # nothing to cancel
                pass
            self._timer_id = self.tree.after(self.delay, self.display_tooltip)

    def display_tooltip(self):
        """Display the tooltip corresponding to the hovered item."""
        item = self.tree.identify_row(self.tree.winfo_pointery() - self.tree.winfo_rooty())
        text = self.tooltip_text.get(item, '')
        self.current_item = item
        if text:
            self.tooltip.configure(text=text)
            self.tooltip.deiconify()
            x = self.tree.winfo_pointerx() + 14
            y = self.tree.winfo_rooty() + self.tree.bbox(item)[1] + self.tree.bbox(item)[3]
            self.tooltip.geometry('+%i+%i' % (x, y))
