#######################################################################
#            Based on RedHat
#######################################################################
  CURRENT_DIR:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
  include $(CURRENT_DIR)/make_redhat.inc

############### Modifications should go below here#####################
