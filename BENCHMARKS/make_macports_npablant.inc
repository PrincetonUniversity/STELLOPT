

#######################################################################
#            Define User Specific Output Paths.
#######################################################################
  # Define the output directory.  This should be an absolute path.
  STELLOPT_HOME = /u/npablant/code/library/stellopt/stellopt_git/bin


#######################################################################
#            Call the standard makefile include for macports.
#######################################################################
  CURRENT_DIR:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
  include $(CURRENT_DIR)/make_macports.inc
