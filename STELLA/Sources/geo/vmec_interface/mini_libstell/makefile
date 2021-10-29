# Typically this makefile is only called from make in the parent directory,
# in which the variables FC and EXTRA_COMPILE_FLAGS are set.

LIBSTELL_OBJ_FILES = ezcdf.o ezcdf_GenPut.o ezcdf_inqvar.o handle_err.o stel_constants.o v3_utilities.o \
	ezcdf_GenGet.o ezcdf_attrib.o ezcdf_opncls.o read_wout_mod.o safe_open_mod.o stel_kinds.o vsvd0.o \
	vparams.o mgrid_mod.o vmec_input.o system_mod.o parse_extension.o getcarg.o vmec_getenv.o

LIBSTELL_TARGET = mini_libstell.a

.PHONY: all clean

all: $(LIBSTELL_TARGET)

%.o: %.f90
	$(FC) $(EXTRA_COMPILE_FLAGS) -DNETCDF -c $<

%.o: %.f
	$(FC) $(EXTRA_COMPILE_FLAGS) -DNETCDF -c $<

%.o: %.F
	$(FC) $(EXTRA_COMPILE_FLAGS) -DNETCDF -c $<

$(LIBSTELL_TARGET): $(LIBSTELL_OBJ_FILES)
	ar rcs $(LIBSTELL_TARGET) $(LIBSTELL_OBJ_FILES)

clean:
	rm -f *.o *.mod *.MOD *~ $(LIBSTELL_TARGET)

# Dependencies:
ezcdf.o: ezcdf_GenPut.o ezcdf_GenGet.o ezcdf_attrib.o ezcdf_opncls.o
ezcdf_GenPut.o: ezcdf_opncls.o ezcdf_inqvar.o
ezcdf_inqvar.o:
handle_err.o:
stel_constants.o: stel_kinds.o
v3_utilities.o: stel_kinds.o
ezcdf_GenGet.o: ezcdf_opncls.o ezcdf_inqvar.o
ezcdf_attrib.o: ezcdf_inqvar.o
ezcdf_opncls.o:
read_wout_mod.o: stel_kinds.o safe_open_mod.o stel_constants.o v3_utilities.o vsvd0.o ezcdf.o mgrid_mod.o vmec_input.o
safe_open_mod.o:
stel_kinds.o:
vsvd0.o: stel_kinds.o
mgrid_mod.o: stel_kinds.o vmec_input.o vsvd0.o system_mod.o safe_open_mod.o ezcdf.o
vparams.o: stel_kinds.o stel_constants.o
vmec_input.o: vparams.o vsvd0.o
system_mod.o:
parse_extension.o:
getcarg.o:
vmec_getenv.o: