# -*- coding: utf-8 -*-
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


Path bar button class
"""


from .constants import add_trace, remove_trace, ttk


class PathButton(ttk.Button):
    """Toggle button class to make the path bar."""

    def __init__(self, parent, variable, value, **kwargs):
        """
        Create a PathButton.

        Like Radiobuttons, only one PathButton in the group (all PathButtons
        sharing the same control variable) can be selected.

        Options:
            * parent: parent widget
            * variable: control variable that the PathButton shares with the
                        other PathButtons in the group (like for Radiobuttons)
            * value: when the PathButton is clicked, the control variable is set
                     to value
            * all ttk.Button options
        """
        kwargs["style"] = "path.tkfilebrowser.TButton"
        kwargs.setdefault("text", "")
        txt = kwargs['text']
        kwargs.setdefault("width", len(txt) + 1 + txt.count('m') + txt.count('M'))
        ttk.Button.__init__(self, parent, **kwargs)
        self.variable = variable
        self.value = value
        self._trace = add_trace(self.variable, "write", self.var_change)
        self.bind("<Button-1>", self.on_press)

    def on_press(self, event):
        """Change the control variable value when the button is pressed."""
        self.variable.set(self.value)

    def get_value(self):
        """Return value."""
        return self.value

    def destroy(self):
        """Remove trace from variable and destroy widget."""
        remove_trace(self.variable, "write", self._trace)
        ttk.Button.destroy(self)

    def var_change(self, *args):
        """Change the state of the button when the control variable changes."""
        self.master.update()
        self.master.update_idletasks()
        if self.variable.get() == self.value:
            self.state(("selected",))
        else:
            self.state(("!selected",))
