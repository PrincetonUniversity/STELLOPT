# Shamelessly stolen from ASCOT4
MACHINE ?= $(shell uname -n | \
	awk '/lazersons-air/{print "macports"}\
	     /Lazersons-MacBook-Air/{print "macports"}\
	     /cori/{print "cori"}\
	     /dawson/{print "pppl_gcc"}\
	     /sunfire/{print "pppl_gcc"}\
	     /kruskal/{print "pppl_gcc"}\
	     /landreman/{print "macports_landreman"}\
	     /Matts-MBP/{print "macports_landreman"}\
	     /Matts-MacBook-Pro/{print "macports_landreman"}\
	     /npablant/{print "macports_npablant"}\
	     /vulcan/{print "vulcan"}\
	     /cobra/{print "cobra"}\
	     /draco/{print "draco"}')

ifndef STELLOPT_PATH
  $(error The environment variable STELLOPT_PATH needs to be set.)
endif

ifneq ($(strip $(MACHINE)),)
	include $(STELLOPT_PATH)/SHARE/make_$(MACHINE).inc
else
  $(error MACHINE was not auto-detected, so you must set the MACHINE environment variable manually.)
endif

# Now include the ALL FILE
include $(STELLOPT_PATH)/SHARE/make_all.inc


