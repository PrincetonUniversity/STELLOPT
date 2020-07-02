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
	mkdir -p $(MYHOME)
  endif



test_make:
	@echo MACHINE is $(MACHINE)
	@echo STELLOPT_HOME is $(STELLOPT_HOME)
	@echo STELLOPT_PATH is $(STELLOPT_PATH)
#	@echo VMEC_DIR is $(VMEC_DIR)
#	@echo BEAMS3D_DIR is $(BEAMS3D_DIR)
#	@echo BOOTSJ_DIR is $(BOOTSJ_DIR)
#	@echo BNORM_DIR is $(BNORM_DIR)
#	@echo BOOZ_DIR is $(BOOZ_DIR)
#	@echo COBRA_DIR is $(COBRA_DIR)
#	@echo DIANGO_DIR is $(DIANGO_DIR)
#	@echo FIELDLINES_DIR is $(FIELDLINES_DIR)
#	@echo JINV_DIR is $(JINV_DIR)
#	@echo MGRID_DIR is $(MGRID_DIR)
#	@echo DKES_DIR is $(DKES_DIR)
#	@echo NEO_DIR is $(NEO_DIR)
#	@echo GENE_DIR is $(GENE_DIR)
#	@echo REGCOIL_DIR is $(REGCOIL_DIR)
#	@echo SFINCS_DIR is $(SFINCS_DIR)
#	@echo MANGO_DIR is $(MANGO_DIR)
	@echo Compiler flags are $(LIBS)
	@echo LIB_LINK is $(LIB_LINK)
