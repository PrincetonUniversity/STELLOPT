# INSTALLATION
This version of STELLOPT makes use of the `make.inc` file to define
compiler options and paths (in analogy to the old setup file).  These
files are stored in the [SHARE](SHARE) folder.  Check to see if one exists
for your computer system.  If it does, just create a symbolic link to
it named `make.inc` in the directory above [SHARE](SHARE).  Otherwise copy 
the `make_pppl.inc` file under a new name. By default, `make_pppl.inc`
assumes that you want to compile all the codes. Please modify your own
copy if needed.

Creating a `make.inc` file is now mandatory, and please DO NOT track 
your personal `make.inc` file via git.

     cp SHARE/make_pppl.inc SHARE/make_pppl_myown.inc
     ln -sf SHARE/make_pppl_myown.inc make.inc

Once a working copy is developed for your computer, you are welcome 
to push it back to the main repository for other users.
Before beginning, it is recommened to create an install branch. To do 
this issue the following commands from main directory:

     git checkout -b install

The -b option tells us to create a new branch (you could also issue
two commands).  Now you'll be on your install branch feel free to
edit files as necessary for your installation.  To build the code 
use the build_all script which will systematically build and compile 
all the files.

For more information, please view [STELLOPT compilation](https://princetonuniversity.github.io/STELLOPT/STELLOPT%20Compilation).


# EDITING CODE
To edit the code please first checkout a new branch with the feature
you'd like.  Then when you want your changes merged into master,
please push your branch to the remote:

     git push origin <branchname>

Then submit a pull request through 
[GitHub](https://github.com/PrincetonUniversity/STELLOPT/pulls).


# PYTHON INTERFACE
A Python interface using CTYPES is also included but it requires a
static shared build of LIBSTELL.  This is still a highly experimental
option.


# REPORT BUGS
If you find any bugs or have any suggestions, please submit an issue 
through [GitHub issues](https://github.com/PrincetonUniversity/STELLOPT/issues).
You are suggested to provide detailed description and assign someone 
when creating an issue.


# DOCUMENTATION
You can find some wiki pages at 
[GitHub pages](https://princetonuniversity.github.io/STELLOPT/).
This is automatically generated from the `markdown` source files 
in the [gh-pages](https://github.com/PrincetonUniversity/STELLOPT/tree/gh-pages) branch. 
Everyone is welcome to contribute to the wiki pages.
