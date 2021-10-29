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


The icons are modified versions of icons from the elementary project
(the xfce fork to be precise https://github.com/shimmerproject/elementary-xfce)
Copyright 2007-2013 elementary LLC.


Constants and functions
"""
import locale
from babel.numbers import format_number
from babel.dates import format_date, format_datetime
from datetime import datetime
import os
from math import log, floor

try:
    import tkinter as tk
    from tkinter import ttk
    from tkinter.messagebox import askyesnocancel, showerror
    from urllib.parse import unquote
except ImportError:
    import Tkinter as tk
    import ttk
    from tkMessageBox import askyesnocancel, showerror
    from urllib import unquote
    import sys
    reload(sys)
    sys.setdefaultencoding('utf8')

PATH = os.path.dirname(__file__)

LOCAL_PATH = os.path.join(os.path.expanduser('~'), '.config', 'tkfilebrowser')

if not os.path.exists(LOCAL_PATH):
    try:
        if not os.path.exists(os.path.join(os.path.expanduser('~'), '.config')):
            os.mkdir(os.path.join(os.path.expanduser('~'), '.config'))
        os.mkdir(LOCAL_PATH)
    except Exception:
        # avoid raising error if the path is not writtable
        pass

RECENT_FILES = os.path.join(LOCAL_PATH, 'recent_files')

# ---  images
if tk.TkVersion < 8.6:
    from PIL.ImageTk import PhotoImage
else:
    PhotoImage = tk.PhotoImage

IM_HOME = os.path.join(PATH, "images", "home.png")
IM_DESKTOP = os.path.join(PATH, "images", "desktop.png")
IM_FOLDER = os.path.join(PATH, "images", "folder.png")
IM_FOLDER_LINK = os.path.join(PATH, "images", "folder_link.png")
IM_NEW = os.path.join(PATH, "images", "new_folder.png")
IM_FILE = os.path.join(PATH, "images", "file.png")
IM_FILE_LINK = os.path.join(PATH, "images", "file_link.png")
IM_LINK_BROKEN = os.path.join(PATH, "images", "link_broken.png")
IM_DRIVE = os.path.join(PATH, "images", "drive.png")
IM_RECENT = os.path.join(PATH, "images", "recent.png")
IM_RECENT_24 = os.path.join(PATH, "images", "recent_24.png")

# ---  translation
try:
    LANG = locale.getdefaultlocale()[0]
except ValueError:
    LANG = 'en'

EN = {}
FR = {"B": "octets", "MB": "Mo", "kB": "ko", "GB": "Go", "TB": "To",
      "Name: ": "Nom : ", "Folder: ": "Dossier : ", "Size": "Taille",
      "Name": "Nom", "Modified": "Modifié", "Save": "Enregistrer",
      "Open": "Ouvrir", "Cancel": "Annuler", "Location": "Emplacement",
      "Today": "Aujourd'hui", "Confirmation": "Confirmation",
      "Error": "Erreur",
      "The file {file} already exists, do you want to replace it?": "Le fichier {file} existe déjà, voulez-vous le remplacer ?",
      "Shortcuts": "Raccourcis", "Save As": "Enregistrer sous",
      "Recent": "Récents", "Recently used": "Récemment utilisés"}
LANGUAGES = {"fr": FR, "en": EN}
if LANG[:2] == "fr":
    TR = LANGUAGES["fr"]
else:
    TR = LANGUAGES["en"]


def _(text):
    """ translation function """
    return TR.get(text, text)


fromtimestamp = datetime.fromtimestamp


def locale_date(date=None):
    return format_date(date, 'short', locale=LANG)


def locale_datetime(date=None):
    return format_datetime(date, 'EEEE HH:mm', locale=LANG)


def locale_number(nb):
    return format_number(nb, locale=LANG)


SIZES = [_("B"), _("kB"), _("MB"), _("GB"), _("TB")]

# ---  locale settings for dates
TODAY = locale_date()
YEAR = datetime.now().year
DAY = int(format_date(None, 'D', locale=LANG))


# ---  functions
def add_trace(variable, mode, callback):
    """
    Add trace to variable.

    Ensure compatibility with old and new trace method.
    mode: "read", "write", "unset" (new syntax)
    """
    try:
        return variable.trace_add(mode, callback)
    except AttributeError:
        # fallback to old method
        return variable.trace(mode[0], callback)


def remove_trace(variable, mode, cbname):
    """
    Remove trace from variable.

    Ensure compatibility with old and new trace method.
    mode: "read", "write", "unset" (new syntax)
    """
    try:
        variable.trace_remove(mode, cbname)
    except AttributeError:
        # fallback to old method
        variable.trace_vdelete(mode[0], cbname)


def get_modification_date(file):
    """Return the modification date of file."""
    try:
        tps = fromtimestamp(os.path.getmtime(file))
    except OSError:
        tps = TODAY
    date = locale_date(tps)
    if date == TODAY:
        date = _("Today") + tps.strftime(" %H:%M")
    elif tps.year == YEAR and (DAY - int(tps.strftime("%j"))) < 7:
        date = locale_datetime(tps)
    return date


def display_modification_date(mtime):
    """Return the modDification date of file."""
    if isinstance(mtime, str):
        return mtime
    tps = fromtimestamp(mtime)
    date = locale_date(tps)
    if date == TODAY:
        date = _("Today") + tps.strftime(" %H:%M")
    elif tps.year == YEAR and (DAY - int(tps.strftime("%j"))) < 7:
        date = locale_datetime(tps)
    return date


def display_size(size_o):
    """Return the size of file."""
    if isinstance(size_o, str):
        return size_o
    if size_o > 0:
        m = int(floor(log(size_o) / log(1024)))
        if m < len(SIZES):
            unit = SIZES[m]
            s = size_o / (1024 ** m)
        else:
            unit = SIZES[-1]
            s = size_o / (1024**(len(SIZES) - 1))
        size = "%s %s" % (locale_number("%.1f" % s), unit)
    else:
        size = "0 " + _("B")
    return size


def key_sort_files(file):
    return file.is_file(), file.name.lower()
