# -*- coding: utf-8 -*-
"""
tkfilebrowser - Alternative to filedialog for Tkinter
Copyright 2017-2018 Juliette Monsel <j_4321@protonmail.com>

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


Functions
"""


from .constants import _
from .filebrowser import FileBrowser


def askopendirname(parent=None, title=_("Open"), **kwargs):
    """
    Return :obj:`''` or the absolute path of the chosen directory.

    Arguments:
    
        parent : Tk or Toplevel instance
            parent window

        title : str
            the title of the filebrowser window

        initialdir : str
            directory whose content is initially displayed

        initialfile : str
            initially selected item (just the name, not the full path)

        filetypes : list :obj:`[("name", "*.ext1|*.ext2|.."), ...]`
          only the files of given filetype will be displayed,
          e.g. to allow the user to switch between displaying only PNG or JPG
          pictures or dispalying all files:
          :obj:`filtypes=[("Pictures", "\*.png|\*.PNG|\*.jpg|\*.JPG'), ("All files", "\*")]`

        okbuttontext : str
            text displayed on the validate button, default is "Open".

        cancelbuttontext : str
            text displayed on the button that cancels the selection, default is "Cancel".

        foldercreation : bool
            enable the user to create new folders if True (default)
    """
    dialog = FileBrowser(parent, mode="opendir", multiple_selection=False,
                         title=title, **kwargs)
    dialog.wait_window(dialog)
    return dialog.get_result()


def askopendirnames(parent=None, title=_("Open"), **kwargs):
    """
    Return :obj:`()` or the tuple of the absolute paths of the chosen directories

    Arguments:
    
        parent : Tk or Toplevel instance
            parent window

        title : str
            the title of the filebrowser window

        initialdir : str
            directory whose content is initially displayed

        initialfile : str
            initially selected item (just the name, not the full path)

        filetypes : list :obj:`[("name", "*.ext1|*.ext2|.."), ...]`
          only the files of given filetype will be displayed,
          e.g. to allow the user to switch between displaying only PNG or JPG
          pictures or dispalying all files:
          :obj:`filtypes=[("Pictures", "\*.png|\*.PNG|\*.jpg|\*.JPG'), ("All files", "\*")]`

        okbuttontext : str
            text displayed on the validate button, default is "Open".

        cancelbuttontext : str
            text displayed on the button that cancels the selection, default is "Cancel".

        foldercreation : bool
            enable the user to create new folders if True (default)
    """
    dialog = FileBrowser(parent, mode="opendir", multiple_selection=True,
                         title=title, **kwargs)
    dialog.wait_window(dialog)
    res = dialog.get_result()
    if not res:  # type consistency: always return a tuple
        res = ()
    return res


def askopenfilename(parent=None, title=_("Open"), **kwargs):
    """
    Return :obj:`''` or the absolute path of the chosen file

    Arguments:
    
        parent : Tk or Toplevel instance
            parent window

        title : str
            the title of the filebrowser window

        initialdir : str
            directory whose content is initially displayed

        initialfile : str
            initially selected item (just the name, not the full path)

        filetypes : list :obj:`[("name", "*.ext1|*.ext2|.."), ...]`
          only the files of given filetype will be displayed,
          e.g. to allow the user to switch between displaying only PNG or JPG
          pictures or dispalying all files:
          :obj:`filtypes=[("Pictures", "\*.png|\*.PNG|\*.jpg|\*.JPG'), ("All files", "\*")]`

        okbuttontext : str
            text displayed on the validate button, default is "Open".

        cancelbuttontext : str
            text displayed on the button that cancels the selection, default is "Cancel".

        foldercreation : bool
            enable the user to create new folders if True (default)
    """
    dialog = FileBrowser(parent, mode="openfile", multiple_selection=False,
                         title=title, **kwargs)
    dialog.wait_window(dialog)
    return dialog.get_result()


def askopenfilenames(parent=None, title=_("Open"), **kwargs):
    """
    Return :obj:`()` or the tuple of the absolute paths of the chosen files

    Arguments:
    
        parent : Tk or Toplevel instance
            parent window

        title : str
            the title of the filebrowser window

        initialdir : str
            directory whose content is initially displayed

        initialfile : str
            initially selected item (just the name, not the full path)

        filetypes : list :obj:`[("name", "*.ext1|*.ext2|.."), ...]`
          only the files of given filetype will be displayed,
          e.g. to allow the user to switch between displaying only PNG or JPG
          pictures or dispalying all files:
          :obj:`filtypes=[("Pictures", "\*.png|\*.PNG|\*.jpg|\*.JPG'), ("All files", "\*")]`

        okbuttontext : str
            text displayed on the validate button, default is "Open".

        cancelbuttontext : str
            text displayed on the button that cancels the selection, default is "Cancel".

        foldercreation : bool
            enable the user to create new folders if True (default)
    """
    dialog = FileBrowser(parent, mode="openfile", multiple_selection=True,
                         title=title, **kwargs)
    dialog.wait_window(dialog)
    res = dialog.get_result()
    if not res:  # type consistency: always return a tuple
        res = ()
    return res


def asksaveasfilename(parent=None, title=_("Save As"), **kwargs):
    """
    Return :obj:`''` or the chosen absolute path (the file might not exist)
    
    Arguments:
    
        parent : Tk or Toplevel instance
            parent window

        title : str
            the title of the filebrowser window

        initialdir : str
            directory whose content is initially displayed

        initialfile : str
            initially selected item (just the name, not the full path)
            
        defaultext : str (e.g. '.png')
            extension added to filename if none is given (default is none)  

        filetypes : list :obj:`[("name", "*.ext1|*.ext2|.."), ...]`
          only the files of given filetype will be displayed,
          e.g. to allow the user to switch between displaying only PNG or JPG
          pictures or dispalying all files:
          :obj:`filtypes=[("Pictures", "\*.png|\*.PNG|\*.jpg|\*.JPG'), ("All files", "\*")]`

        okbuttontext : str
            text displayed on the validate button, default is "Open".

        cancelbuttontext : str
            text displayed on the button that cancels the selection, default is "Cancel".

        foldercreation : bool
            enable the user to create new folders if True (default)
    """
    dialog = FileBrowser(parent, mode="save", title=title, **kwargs)
    dialog.wait_window(dialog)
    return dialog.get_result()
