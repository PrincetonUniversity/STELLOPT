# This top level makefile is only used to setup the build directories.
# In general this should not be called directly, instead use the
# build_all script.
#
#  This top level makefile includes
#      make.inc includes
#           make_$(MACHINE).inc (from SHARE)  ... no external makefiles
# and       make_all.inc        (from SHARE)
#
#

include make.inc

FLAG_CALLED_FROM_BUILD_ALL ?= false

.PHONY: release debug clean_release clean_debug static_release shared_release

release: pre_build

debug: pre_build

clean_release: pre_build

clean_debug: pre_build

static_release: pre_build

shared_release: pre_build

pre_build:
ifeq "$(FLAG_CALLED_FROM_BUILD_ALL)" "false"
	@echo 'To compile STELLOPT use the build_all script instead of calling make directly.'
	@(exit 1)
endif
ifeq ($(wildcard $(MYHOME)),)
	@echo 'Creating STELLOPT output directory.'
	mkdir -p $(MYHOME)
endif

pystel: libstell$(SHARED_EXT)
	@echo 'Building pySTEL'
	@cd pySTEL; python3 setup.py install --user

libstell$(SHARED_EXT):
	@cd LIBSTELL; make shared_release

test_make:
	@echo 'Directories and flags for build.'
	@echo '--------------------------------'
	@echo MACHINE is $(MACHINE)
	@echo STELLOPT_HOME is $(STELLOPT_HOME)
	@echo STELLOPT_PATH is $(STELLOPT_PATH)
	@echo '--------------------------------'
	@echo " LAEOPT:          " $(LAEOPT) " DIR " $(AEOPT_DIR)
	@echo " LCOILOPT:        " $(LCOILOPT) " DIR "  $(COILOPTPP_DIR)
	@echo " LGENE:           " $(LGENE) " DIR " $(GENE_DIR)
	@echo " LMANGO:          " $(LMANGO) " DIR " $(MANGO_DIR)
	@echo " LREGCOIL:        " $(LREGCOIL) " DIR " $(REGCOIL_DIR)
	@echo " LSINFCS:         " $(LSFINCS) " DIR " $(SFINCS_DIR)
	@echo " LTERPSICHORE:    " $(LTERPSICHORE) " DIR " $(TERPSICHORE_DIR)
	@echo " LTRAVIS:         " $(LTRAVIS) " DIR " $(TRAVIS_DIR)
	@echo '--------------------------------'
	@echo Precompiler flags are $(PRECOMP)
	@echo '--------------------------------'
	@echo Compiler $(COMPILE)
	@echo Compiler Free $(COMPILE_FREE)
	@echo Compiler modules are $(MOD1_PATH)
	@echo Compiler release flags are $(FLAGS_R)
	@echo Compiler debug flags are $(FLAGS_D)
	@echo '--------------------------------'
	@echo Linker is $(LINK)
	@echo Linker flags are $(LIB_LINK)
	@echo '--------------------------------'
	@echo `python --version` `which python`
	@echo '--------------------------------'
