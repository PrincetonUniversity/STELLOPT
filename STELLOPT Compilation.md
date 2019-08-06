STELLOPT Suite Compilation
==========================

------------------------------------------------------------------------

Overview
--------

The STELLOPT Suite compiles a subset of included codes and can link to
precompiled codes (expanding it\'s capabilities). STELLOPT is maintained
under Git. The repository is maintained at Bitbucket.

Installation Workflow
---------------------

In general users must make slight modifications to the files in the repo
before compiling. For this reason we ask that users create an
\'install\' branch based on the branch from which they wish to work. The
**master** branch contains the cutting edge version of the code while
**vXXX** are the release versions of the code. The general workflow
should be:

1.  Clone Repository
2.  Checkout an \'install\' branch locally
3.  Repoint make.inc to the proper SHARE/make\_XXX.inc file (and make
    one if necessary).
4.  Issue the command make from the root directory.

Development Workflow
--------------------

For those wishing to develop the code please create a new branch based
on **master**. Before pushing this branch back to bitbucket, do a pull
and merge of the master branch into your branch to verify that updates
and changes work. Then push you branch to the Bitbucket server. Once you
done this, issue a pull request to have your branch merged with master.
All merges into master should be tracked through the Bitbucket interface
(save small things which do not directly affect compiling).

Specific Compilation Issues
---------------------------

[Compiling on the PPPL Cluster](STELLOPT Compilation at PPPL)

[Compiling on an Apple Macintosh under OSX](STELLOPT Compilation OSX)

[Compiling at NERSC](STELLOPT Compilation at NERSC)

[Compiling on a Redhat/CentOS (RPM) Machine](STELLOPT Compilation CentOS)

[Compiling on the RZG-MPE Cluster (Hydra, Draco, Cobra)](STELLOPT Compilation on Hydra (RZG-MPG))

[Compiling on an Ubuntu (Debian) Machine](STELLOPT Compilation Ubuntu)
