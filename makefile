# This top level makefile is only used to setup the build directories.
# In general this should not be called directly, instead use the
# build_all script.

include make.inc

.PHONY: release debug clean_release clean_debug static_release shared_release

release: pre_build

debug: pre_build

clean_release: pre_build

clean_debug: pre_build

static_release: pre_build

shared_release: pre_build

pre_build:
  ifneq ($(FLAG_CALLED_FROM_BUILD_ALL),true)
	@echo 'To compile STELLOPT use the build_all script instead of calling make directly.'
	@(exit 1)
  endif
  ifeq ($(wildcard $(MYHOME)),)
	@echo 'Creating STELLOPT output directory.'
	mkdir $(MYHOME)
  endif



test_make:
	@echo MACHINE is $(MACHINE)
	@echo STELLOPT_HOME is $(STELLOPT_HOME)
	@echo FC_LINKER_FLAGS is $(FC_LINKER_FLAGS)
	@echo FFLAGS is $(FFLAGS)
	@echo FLINKER is $(FLINKER)
	@echo PETSC_LIB is $(PETSC_LIB)
	@echo SFINCS_LIB is $(SFINCS_LIB)
