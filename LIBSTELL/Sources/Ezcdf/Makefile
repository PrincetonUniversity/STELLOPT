#Makefile to build libezcdf.a
#

ifneq ($(MAKELEVEL),0)
# if OBJ was defined for main Makefile one level down
ifeq ("${OBJ}",".")
	override OBJ=..
endif
endif

SHR = ../share/Make.flags
ifneq ("$(wildcard $(CODESYSDIR)/source/misc/makeflags.mk)","")
	SHR = $(CODESYSDIR)/source/misc/makeflags.mk
endif
-include ../share/Make.local 
include $(SHR)


# for HP : Optimizer does not handle f90 interfaces
ifeq ($(MACHINE),HP)
	FFLAGS = -c
endif
# also problem on AIX
ifeq ($(MACHINE),RS6000)
	FFLAGS = -c
endif

ifdef DEBUG
	FFLAGS = $(DFFLAGS)
	CFLAGS = $(DCFLAGS)
endif
 
OBJDIR = $(OBJ)/obj/ezcdf
MDIR = $(OBJ)/mod

srcdir = $(shell pwd)
 
# for TRANSP: use <name>.a instead of unix convention: lib<name>.a
#
ifdef CODESYSDIR
	ARC = $(OBJ)/lib/ezcdf.a
	MODFLAGS = $(MFLAGS)$(OBJ)/mod
else
	ARC = $(OBJ)/lib/libezcdf.a
endif

# don't rebuild library if using installed one in $PREFIX 
ifeq ($(MAKELEVEL),0)
	THISLIB=$(ARC)
endif
 
# f90 modules
 
Mobjs = ezcdf_inqint.o ezcdf_genput.o ezcdf_genget.o ezcdf.o
MODS = $(foreach file,$(Mobjs),$(ARC)($(file)))


ifeq ($(MACHINE),CRAY)
   INCFLAGS := $(INCFLAGS) -I$(NETCDF_DIR)/include
   CINCL := -I$(NETCDF_DIR)/include
else
   INCFLAGS := $(INCFLAGS) $(IOPT)$(LIBROOT)/include
endif

ifeq ($(MODUP),Y)
 MODS0=$(foreach m,$(Mobjs),$(shell  echo $(m) | tr 'a-z' 'A-Z'))
 MODULES=$(foreach m,$(MODS0),$(subst .O,.$(MODEXT),$(m)))
else
 MODULES = $(foreach m,$(Mobjs),$(subst .o,.$(MODEXT),$(m)))
endif 


# library members
MEM0 = $(filter-out eztest.f90,$(filter-out ezcdf%,$(wildcard *.f90))) ezcdf_open.f90 ezcdf_close.f90

MEM = $(foreach m,$(MEM0),$(ARC)($(subst .f90,.o,$(m))))
 
.PHONY: clean realclean test checkdir
 
# make sure Object directory exists
 
libs: libezcdf

libezcdf:  FORCE $(ARC)
 
all:     libezcdf exec
 
$(OLDLIB): timestamp.inf
	@echo "--- DETECTED libezcd.a source update"
	@echo "--- Re-Making $(ARC)"
	@$(MAKE) libs
 
FORCE:
	@test -d $(OBJ)/lib || mkdir -p $(OBJ)/lib
	@test -d $(OBJ)/mod || mkdir -p $(OBJ)/mod
	@test -d $(OBJDIR) || mkdir -p $(OBJDIR)
 
 
$(OBJDIR)/%.o: %.f90
	$(FC90)  $(FFLAGS) $(MODFLAGS) $(INCFLAGS) $(OUT) $@ $<

$(OBJDIR)/%.o: %.cc
	$(CXX) $(CXXFLAGS) $(DEFS) $(CINCL) -o $@ $<
 
$(ARC)(%.o): %.f90
	$(FC90)  $(FFLAGS) $(MODFLAGS) $(INCFLAGS) $(OUT) $(OBJDIR)/$*.o $<
	$(mv_mod)
	$(make_ar)
 
 
$(ARC): $(MODS) $(MEM)
	@$(RANLIB) $@ > /dev/null
 
exec: ckdirs $(OBJ)/test/eztest $(OBJ)/test/eztest2
 
ckdirs: FORCE
	@test -d $(OBJ)/test || mkdir -p $(OBJ)/test
	@(cd $(OBJ)/test; test -f EZsample.cdf || ln -s $(srcdir)/EZsample.cdf EZsample.cdf )
 
$(OBJ)/test/eztest: $(OBJDIR)/eztest.o $(THISLIB) 
	$(FC90) $(LDFLAGS) -o $@ $< -L$(OBJ)/lib $(LLOC) -lezcdf -L$(USRLIB) $(NETCDF)
 
$(OBJ)/test/eztest2: $(OBJDIR)/eztest2.o $(THISLIB)
	$(CXX) $(LDFLAGS) -o $@ $< -L$(OBJ)/lib $(LLOC) -lezcdf -L$(USRLIB) $(NETCDF) $(XLIBS)
 
install:
	@test -d $(LIBDIR) || mkdir -p $(LIBDIR)
	@test -d $(MODDIR) || mkdir -p $(MODDIR)
	@if test -f $(ARC); then \
	  umask 133; cp $(ARC) $(LIBDIR)/; \
	else test -f $(LIBDIR)/libezcdf.a || echo ">>> can not find $(ARC)";fi
	@(cd $(MDIR); umask 133; cp $(MODULES) $(MODDIR)/)
	@test -d $(MANDIR)/man3 || mkdir -p $(MANDIR)/man3
	@umask 133; cp  ezcdf.3 $(MANDIR)/man3/

uninstall:
	rm -f $(LIBDIR)/libezcdf.a
	rm $(MANDIR)/man3/ezcdf.3
	(cd $(MODDIR); rm -f $(MODULES)) 
clean:
	@rm -f $(OBJDIR)/*
	@rm -f $(OBJ)/test/eztest* 
 
realclean: clean
	@rm -f $(ARC)
	@(cd $(MDIR); rm -f $(MODULES))
