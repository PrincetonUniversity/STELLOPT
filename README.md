# INSTALLATION
This version of STELLOPT makes use of `make_***.inc` files in the `SHARE`
subdirectory to set various options at compile time.  The `make.inc` file
in the main directory is a script which picks these files based on what
is returned by the `uname -n` command on your machine.  This can be
overridden by setting the environment variable `MACHINE` to the file you
wish to use.  For example to use `SHARE/make_bobdole.inc` you would set
`MACHINE=bobdole` before calling the build scripts.  It is also
important to set `STELLOPT_PATH` to the path to your current directory
where you've pulled the repository.

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
