=======================================================================
=====     How to compile STELLOPT                                   ===
=======================================================================

This version of STELLOPT makes use of the make.inc files to define
compiler options and paths (in analogy to the old setup file).  These
files are stored in the SHARE folder.  Check to see if one exists
for your computer system.  If it does just create a symbolic link to
it named make.inc in the directory above SHARE.  Otherwise
copy the make_pppl.inc file under a new name and modify as needed.
Once a working copy is developed for your computer make sure to push
it back to the main repository for other users.

To build the code use the build_all script which will systematically
build and compile all the files.

A Python interface using CTYPES is also included but it requires a
static shared build of LIBSTELL.  This is still a highly experimental
option.
