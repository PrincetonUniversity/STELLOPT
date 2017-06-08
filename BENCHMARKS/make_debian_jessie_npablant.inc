
#######################################################################
#            Find the source paths of the STELLOPT distribution.
#######################################################################
  CURRENT_DIR:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
  STELLOPT_DIR:=$(abspath $(CURRENT_DIR)/../)


#######################################################################
#            Define User Specific Output Paths.
#######################################################################
  # Define the output directory.  This should be an absolute path.
  STELLOPT_HOME = $(STELLOPT_DIR)/bin


#######################################################################
#            Call the standard makefile include for macports.
#######################################################################
  include $(CURRENT_DIR)/make_debian_jessie.inc
