#Type "make release(debug)" at command line to make release(debug) version  
#DO NOT alter the directory structure or move this makefile (must be in first subdirectory below top_level_dir)
  include ../make.inc
  DBG = Debug
  REL = Release
  REL_PATH = $(PWD1)/$(REL)
  DBG_PATH = $(PWD1)/$(DBG)
  EXE_F = xvmec2000
  MAKE = makevmec2000
  LIB_DIR = $(MYHOME)/libstell_dir

.PHONY: test_vmec

test_vmec:
	@echo "Beginning testing of VMEC2000"
	@cd VMEC_TEST
	@rm -f jxbout* mercier* threed1* wout*
	@srun -n 8 -p dawson -t 20:00:00 ~/bin/xvmec2000 input.DSHAPE
	@./compare_DSHARE.py
	@echo "Tesing of VMEC2000 complete"
	@echo ""

