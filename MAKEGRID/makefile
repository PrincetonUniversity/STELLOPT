#Type "make release(debug)" at command line to make release(debug) version  
#DO NOT alter the directory structure or move this makefile (must be in first subdirectory below top_level_dir)
  include ../make.inc
  DBG = Debug
  REL = Release
  REL_PATH = $(PWD1)/$(REL)
  DBG_PATH = $(PWD1)/$(DBG)
  EXE_F = xgrid
  MAKE = makemgrid
  LIB_DIR = $(MYHOME)/libstell_dir

.PHONY: release debug clean_release clean_debug static_release

release:
	@echo "Beginning release compilation and linking of makegrid"
	@cd $(REL); rm -f *.f *.f90 *.c; rm -f $(MAKE); ln -fs ../$(MAKE) .; make -f $(MAKE) $(EXE_F) FLAGS="$(FLAGS_R) $(MOD1_PATH)" TYPE="release" LOCTYPE="Release"
	@rm -f $(REL)\*.f; rm -f $(REL)\*.f90
	@rm -f $(MYHOME)/$(EXE_F)
	@-(umask 022; rm -f $(MYHOME)/$(EXE_F); ln -fs $(REL_PATH)/$(EXE_F) $(MYHOME) 2> /dev/null; umask 077)
	@echo "Release version - xgrid - is now updated"
	@echo "To free up space, delete all *.o files from Release"
	@echo ""

debug:
	@echo "Beginning debug compilation and linking of makegrid"
	@cd $(DBG); rm -f *.f *.f90 *.c; rm -f $(MAKE); ln -fs ../$(MAKE) .; make -f $(MAKE) $(EXE_F) FLAGS="$(FLAGS_D) $(MOD1_PATH)" TYPE="debug" LOCTYPE="Debug" 
	@rm -f $(MYHOME)/$(EXE_F)
	@-(umask 022; rm -f $(MYHOME)/$(EXE_F); ln -fs $(DBG_PATH)/$(EXE_F) $(MYHOME) 2> /dev/null; umask 077)
	@echo "Debug version - xgrid - is now updated"
	@echo "To free up space, delete all *.o files from Debug"
	@echo ""

clean_release: 
	@echo "CLEAN REMAKE IN PROGRESS"
	@cd $(REL); make -f $(MAKE) clean
	@make release

clean_debug: 
	@echo "CLEAN REMAKE IN PROGRESS"
	@cd $(DBG); make -f $(MAKE) clean
	@make debug
      
static_release:      
	@echo "MAKING STATIC (RELEASE) EXECUTABLE"
	@cd $(REL); rm -f *.f *.f90 *.c; rm -f $(MAKE); ln -fs ../$(MAKE) .; make -f $(MAKE) $(EXE_F) FLAGS="$(FLAGS_R) $(MOD1_PATH)" SFLAGS="" TYPE="release"
	@rm -f $(MYHOME)/$(EXE_F)
	@-(umask 022; rm -f $(MYHOME)/$(EXE_F); ln -fs $(REL_PATH)/$(EXE_F) $(MYHOME) 2> /dev/null; umask 077)
	@echo "Static Release version - xgrid - is now updated"
	@echo "To free up space, delete all *.o files from Release"
	@echo ""
