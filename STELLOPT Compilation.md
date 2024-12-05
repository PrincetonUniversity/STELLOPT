STELLOPT Suite Compilation
==========================

------------------------------------------------------------------------

Overview
--------

The STELLOPT Suite compiles a subset of included codes and can link to
precompiled codes (expanding it\'s capabilities). STELLOPT is maintained
using git. The repository is maintained at GitHub.

To use the MANGO library of optimization algorithms in STELLOPT,
you should first obtain and build MANGO following the [instructions here](https://hiddensymmetries.github.io/mango/gettingStarted.html).

Get source code
----------

The STELLOPT code is now fully open-sourced with MIT license. You can view the source code at https://github.com/PrincetonUniversity/STELLOPT.

Installation Workflow
---------------------
The **master** branch contains the cutting edge version of the code while
**develop** has the latest testing capabilities.The general workflow
should be:

1.  Clone Repository
2.  Set an environment variable `STELLOPT_PATH` equal to the path to the repository folder (for example `~/src/STELLOPT/`)
3.  Some cluster and devices require the `MACHINE` variable to be set (specific links below).
4.  Issue the command `./build_all` to build all the codes.

To build a specific subcode enter the directory of the code and issue the command `make clean_release` to rebuild the whole code with optimization flags, `make clean_debug` to build the whole code for debugging, `make release` to only build files which have changed, `make debug` to only build files which have changed with debugging on, or `make` which aliases to `make release`.

Precompiler Flags
-----------------

Certain features of codes are controled via precompiler flags which are not by default turned on.  The following is a list of precompiler flags one can add to their `SHARE/make_machine.inc` file by adding lines such as `PRECOMP+= -DHDF5_PAR`

- `-DNETCDF`:   NetCDF output (VMEC2000, MAKEGRID)
- `-DLHDF5`:     HDF5 output (FIELDLINES, BEAMS3D)
- `-DHDF5_PAR`: Parallel HDF5 routines (FIELDLINES, BEAMS3D)
- `-DNAG`:     NAG Algorithms (FIELDLINES, BEAMS3D)
- `-DIMAS`:     IMAS Interface (VMEC2000, BEAMS3D)
- `-D_ANIMEC`:  Anisotropic VMEC (VMEC2000)
- `-D_FLOW`:  Toroidal flow VMEC (VMEC2000)
- `-DCHI_FORCE`:  Alternative VMEC forces (VMEC2000)
- `-DB3D_COLLOP_NRL19IE`:  Ion/Electron NRL19 Collision operator (BEAMS3D)
- `-DB3D_COLLOP_NUBEAM`:  NUBEAM Collision operator (BEAMS3D)
- `-DB3D_VEL_DIFFUSION`:  Velocity diffusion (BEAMS3D)

There are others which are control via the presence of secondary codes, only the above should be manually set and only if you know what you are doing.

Development Workflow
--------------------

For those wishing to develop the code please create a new branch based
on **develop** . Before pushing this branch back to bitbucket, do a pull
and merge of the **develop** branch into your branch to verify that updates
and changes work. Then push you branch to the GitHub server. Once you
done this, issue a pull request to have your branch merged with **develop**.
All merges into **develop** should be tracked through the GitHub interface
(save small things which do not directly affect compiling).

Specific Compilation Issues
---------------------------

*If you find the instructions are out dated, please raise an issue. Direct contributions to the documentation on the **gh-pages** branch are particularly welcomed.*

[Compiling on an Apple Macintosh under OSX](STELLOPT Compilation OSX)

[Compiling on the CINECA Cluster (Marconi)](STELLOPT Compilation CINECA)

[Compiling on an Debian Machine](STELLOPT Compilation Debian)

[Compiling with Docker](STELLOPT Compilation Docker)

[STELLOPT Compilation on the EUROfusion Gateway](STELLOPT Compilation on the EUROfusion Gateway)

[Compiling on the ITER Cluster](STELLOPT Compilation at ITER)

[Compiling on the MPCDF Clusters (Raven, Viper)](STELLOPT Compilation on MPCDF)

[Compiling at NERSC](STELLOPT Compilation at NERSC)

[Compiling on the PPPL Cluster](STELLOPT Compilation at PPPL)

[Compiling on a Redhat Machine](STELLOPT Compilation Redhat)

