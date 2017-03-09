=======================================================================
=====     How to compile STELLOPT                                   ===
=======================================================================

If you can read this file then you've unzipped the stellinstall_xxx.zip
file.  Various cshrc files can be found in the BENCHMARKS directory 
along with test input files.  You will need to 'source' the
apporpriate cshrc before compiling.  Additionally, code based
environment varialbes (GENE,COILOTP++,etc.) need to be properly set
for you system.

The next step will be to run the setup script
./setup
It will prompt you for a directory in which to place symbolic links
to all compiled packages.  The default is ~/bin/, however if you'd
like to have multiple versions to work with you can enter in a
specific name.  You should indicate you'd like to build a
Release verison.  You can should then answer yes to all questions 
and specify that you'd like Clean build.  Note that as the compiler
builds the codes it will skip over any codes which were not included
in your zip file.
