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


The icons are modified versions of icons from the elementary project
(the xfce fork to be precise https://github.com/shimmerproject/elementary-xfce)
Copyright 2007-2013 elementary LLC.


Recent files management
"""


class RecentFiles:
    """Recent files manager."""
    def __init__(self, filename, nbmax=30):
        """
        Create a recent file manager.

        Options:
            * filename: file where the recent file list is read/saved
            * nbmax: maximum number of recent files to remember
        """
        self._filename = filename
        self.nbmax = nbmax
        self._files = []  # most recent files first
        try:
            with open(filename) as file:
                self._files = file.read().splitlines()
        except Exception:
            pass

    def get(self):
        """Return recent file list."""
        return self._files

    def add(self, file):
        """Add file to recent files."""
        if file not in self._files:
            self._files.insert(0, file)
            if len(self._files) > self.nbmax:
                del(self._files[-1])
        else:
            self._files.remove(file)
            self._files.insert(0, file)
        try:
            with open(self._filename, 'w') as file:
                file.write('\n'.join(self._files))
        except Exception:
            # avoid raising errors if location is read-only or invalid path
            pass
