# Shamelessly stolen from ASCOT4
MACHINE ?= $(shell uname -n | \
	awk '/dion-ubuntu/{print "ubuntu"}\
	     /lazersons-air/{print "macports"}\
	     /Lazersons-MacBook-Air/{print "macports"}\
	     /cori/{print "cori"}\
	     /dawson/{print "pppl_gcc"}\
	     /sunfire/{print "pppl_gcc"}\
	     /kruskal/{print "pppl_gcc"}\
	     /landreman/{print "macports_landreman"}\
	     /npablant/{print "macports_npablant"}\
	     /vulcan/{print "vulcan"}\
	     /Matts-MBP/{print "macports_landreman"}\
	     /Matts-MacBook-Pro/{print "macports_landreman"}\
	     /raven/{print "raven"}\
	     /cobra/{print "cobra"}\
	     /draco/{print "draco"}')

ifneq ($(strip $(MACHINE)),)
	include $(STELLOPT_PATH)/SHARE/make_$(MACHINE).inc
endif

# Now include the ALL FILE
include $(STELLOPT_PATH)/SHARE/make_all.inc


