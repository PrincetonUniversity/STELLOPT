# The variable LIBSTELL_DIR should either be "mini_libstell", if you use this reduced version of libstell
# that comes packaged with this repository, or else it should point to a directory containing libstell .mod files
# elsewhere on your system.
LIBSTELL_DIR = mini_libstell
#LIBSTELL_DIR=/Users/mattland/stellopt/LIBSTELL/Release

# The variable LIBSTELL_FOR_STELLA should either be "mini_libstell/mini_libstell.a", if you use this reduced version of libstell
# that comes packaged with this repository, or else it should point to a libstell.a library elsewhere on your system.
LIBSTELL_FOR_STELLA = mini_libstell/mini_libstell.a
#LIBSTELL_FOR_STELLA=/Users/mattland/stellopt/LIBSTELL/Release/libstell.a

# End of system-dependent variable assignments

# If F90FLAGS has not been set (e.g. running makefile from geo/vmec_interface directory), then set it to
# 'empty' so we can check for it later
F90FLAGS ?= empty

# Load up the compiler flags from the main Stella makefile. We need to include the NETCDF paths for mini_libstell 
EXTRA_COMPILE_FLAGS=$(F90FLAGS) $(NETCDF_INC)
EXTRA_LINK_FLAGS=$(NETCDF_LIB)

export

.PHONY: all clean

all: test_vmec_to_stella_geometry_interface


include makefile.depend

%.o: %.f90
	$(FC) $(EXTRA_COMPILE_FLAGS) -I $(LIBSTELL_DIR) -c $<

%.o: %.f
	$(FC) $(EXTRA_COMPILE_FLAGS) -I $(LIBSTELL_DIR) -c $<

ifeq ($(F90FLAGS),empty)
test_vmec_to_stella_geometry_interface: 
	$(MAKE) vmec -C ../..
else
test_vmec_to_stella_geometry_interface: $(OBJ_FILES)
	$(FC) -o test_vmec_to_stella_geometry_interface $(OBJ_FILES) $(LIBSTELL_FOR_STELLA) $(EXTRA_LINK_FLAGS)
endif

ifeq ($(F90FLAGS),empty)
mini_libstell/mini_libstell.a:
	$(MAKE) mini_libstell.a -C ../..
else
mini_libstell/mini_libstell.a:
	$(MAKE) -C mini_libstell
endif

clean:
	rm -f *.o *.mod *.MOD *~ test_vmec_to_stella_geometry_interface
	cd mini_libstell; rm -f *.o *.mod *.MOD *.a
