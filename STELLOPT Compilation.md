STELLOPT Suite Compilation
==========================

------------------------------------------------------------------------

Overview
--------

The STELLOPT Suite compiles a subset of included codes and can link to
precompiled codes (expanding it\'s capabilities). STELLOPT is maintained
using git. The repository is maintained at GitHub.

Get source code
----------

If you want to get access to the source code, please send your GitHub account 
name to Dr. Caoxiang Zhu (czhu[_at_]pppl[_dot_]gov) with some brief justifications.
Once you obtained the access, you can visit the GitHub page at https://github.com/PrincetonUniversity/STELLOPT.

Installation Workflow
---------------------
The **master** branch contains the cutting edge version of the code while
**develop** has the latest testing capabilities.The general workflow
should be:

1.  Clone Repository
2.  Repoint make.inc to the proper SHARE/make\_XXX.inc file (or make
    one if necessary).
3.  Issue the command make from the root directory.

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

[Compiling on the PPPL Cluster](STELLOPT Compilation at PPPL)

[Compiling on an Apple Macintosh under OSX](STELLOPT Compilation OSX)

[Compiling at NERSC](STELLOPT Compilation at NERSC)

[Compiling on a Redhat/CentOS (RPM) Machine](STELLOPT Compilation CentOS)

[Compiling on the RZG-MPE Cluster (Hydra, Draco, Cobra)](STELLOPT Compilation on Hydra (RZG-MPG))

[Compiling on the CINECA Cluster (Marconi)](STELLOPT Compilation CINECA)

[Compiling on an Ubuntu (Debian) Machine](STELLOPT Compilation Ubuntu)
