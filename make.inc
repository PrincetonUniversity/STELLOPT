# Shamelessly stolen from ASCOT4
MACHINE ?= $(shell uname -n | \
	awk '/lazersons-air/{print "macports"}\
	     /cori/{print "cori"}\
	     /dawson/{print "pppl"}\
	     /sunfire/{print "pppl"}\
	     /kruskal/{print "pppl"}\
	     /landreman/{print "macports_landreman"}\
	     /npablant/{print "macports_npablant"}\
	     /vulcan/{print "vulcan"}\
	     /cobra/{print "cobra"}\
	     /draco/{print "draco"}')

ifneq ($(strip $(MACHINE)),)
	include $(STELLOPT_PATH)/SHARE/make_$(MACHINE).inc
endif


