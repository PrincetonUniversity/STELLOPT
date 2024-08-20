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

Once you configure the environment variables, you can begin to compile the code by

     ./build_all 

By default, it will compile all the codes in the folder with the option `clean_release`.
You can customize the script by providing additional commandline arguments, like `./build_all -o release -j 4 STELLOPTV2`. To check the options, type `./build_all -h`.

For more information, please view [STELLOPT compilation](https://princetonuniversity.github.io/STELLOPT/STELLOPT%20Compilation).

# DOCKER IMAGE
We provide a Dockerfile for building a STELLOPT image based on `debian:latest`. This is essentially the same base image we use for running our workflows.  There is also a very old Docker image available at:
https://hub.docker.com/r/zhucaoxiang/stellopt.

# CITING CODE
```
@misc{ doecode_12551,
title = {STELLOPT},
author = {Lazerson, Samuel and Schmitt, John and Zhu, Caoxiang and Breslau, Joshua and STELLOPT Developers, All},
abstractNote = {The STELLOPT code is designed to optimize 3D MHD equilibria to a set of target physics parameters encompassing stellarator design and 3D equilibrium reconstructions.},
url = {https://doi.org/10.11578/dc.20180627.6},
howpublished = {[Computer Software] \url{https://doi.org/10.11578/dc.20180627.6}},
year = {2020},
month = {may}
}
```

# EDITING CODE
Once a working copy is developed for your computer, you are welcome 
to push it back to the main repository for other users.
You can ask for a write permission and push your changes into a remote branch.
Then submit a pull request through 
[GitHub](https://github.com/PrincetonUniversity/STELLOPT/pulls).

Or you can share your changes via [git fork](https://docs.github.com/en/github/getting-started-with-github/fork-a-repo), which doesn't require you have the write permission.

GitHub actions will be triggered by each push and pull-request. The code will be compiled with debug option and you can check the full log at [GitHub actions](https://github.com/PrincetonUniversity/STELLOPT/actions).

# UTILITIES
A Python interface using CTYPES is also included but it requires a
static shared build of LIBSTELL. We essentially attempt to provide a
class for each code, and a separate utility to call for easy file manipulation
and plotting. 2D plotting in handled by Matplotlib while 3D is handled by VTK.
We provide a wrapper class for VTK to ease integration. You can find the python scripts at `./pySTEL`.

There are also some other plotting packages developed by group members.
Here are some examples.

  - [matlabVMEC](https://github.com/lazersos/matlabVMEC): MATLAB packages for various codes.
  - [coilpy](https://github.com/zhucaoxiang/CoilPy): Python package for toroidal Fourier surface, coils, FOCUS and STELLOPT.

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
You are welcome to report/fix/update the wiki pages.
