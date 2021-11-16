#!/bin/sh

# Make shortcuts to path folders
export STELLA='/home/hanne/Dropbox/stella'
export STELLAPY='/home/hanne/Dropbox/stella/stellapy'

# Run the script
/usr/bin/python3 $STELLAPY/GUI/stella_GUI.py >> ./crash.log 2>&1
