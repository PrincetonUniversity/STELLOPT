#Type "make release(debug)" at command line to make release(debug) version  
#DO NOT alter the directory structure or move this makefile (must be in first subdirectory below top_level_dir)
  include ../make.inc
  LIB_DIR = $(MYHOME)/libstell_dir
  BIN_DIR = ~/bin_cori
  VMEC_TEST_DIR = VMEC_TEST
  DIAGNO_TEST_DIR = DIAGNO_TEST
  STELLOPT_TEST_DIR = STELLOPT_TEST

.PHONY: test_vmec

test_vmec:
	@echo "Beginning testing of VMEC2000"
	@cd $(VMEC_TEST_DIR); rm -f jxbout* mercier* threed1* wout*; 
	@echo $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xvmec2000 input.DSHAPE
	@cd $(VMEC_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xvmec2000 input.DSHAPE
	@cd $(VMEC_TEST_DIR) && ./compare_DSHAPE.py
	@echo ""
	@echo $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xvmec2000 input.SOLOVEV
	@cd $(VMEC_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xvmec2000 input.SOLOVEV
	@cd $(VMEC_TEST_DIR); ./compare_SOLOVEV.py
	@echo ""
	@echo $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xvmec2000 input.HELIOTRON
	@cd $(VMEC_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xvmec2000 input.HELIOTRON
	@cd $(VMEC_TEST_DIR); ./compare_HELIOTRON.py
	@echo ""
	@echo $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xvmec2000 input.ATF
	@cd $(VMEC_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xvmec2000 input.ATF
	@cd $(VMEC_TEST_DIR); ./compare_ATF.py
	@echo ""
	@echo $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xvmec2000 input.QAS
	@cd $(VMEC_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xvmec2000 input.QAS
	@cd $(VMEC_TEST_DIR); ./compare_QAS.py
	@echo "Tesing of VMEC2000 complete"
	@echo ""

test_diagno:
	@echo "Beginning testing of DIAGNO"
	@cd $(DIAGNO_TEST_DIR); rm -f diagno_bench.* diagno_flux.* diagno_surf.* test_loops_b.* test_loops_j.*; 
	@echo $(MPI_RUN) $(MPI_RUN_OPTS) ~/bin/xdiagno -vmec bigtok -bench
	@cd $(DIAGNO_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS) ~/bin/xdiagno -bench -vmec bigtok
	@cd $(DIAGNO_TEST_DIR); ./compare_bigtok.py
	@echo ""

test_stellopt:
	@echo "Beginning testing of STELLOPT"
	@cd $(STELLOPT_TEST_DIR)
	@echo $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.IOTA_LMDIF
	@cd $(STELLOPT_TEST_DIR)/IOTA_LMDIF; $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.IOTA_LMDIF
	@echo ""
	@echo $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.TOK_R0_DELTA
	@cd $(STELLOPT_TEST_DIR)/TOK_R0_DELTA; $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.TOK_R0_DELTA
	@echo ""
	@echo $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.TOK_R0_RHO
	@cd $(STELLOPT_TEST_DIR)/TOK_R0_RHO; $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.TOK_R0_RHO
	@echo ""
	@echo $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.vboot_tokamak_bootsj
	@cd $(STELLOPT_TEST_DIR)/vboot_tokamak_bootsj; $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.vboot_tokamak_bootsj
	@if [ $LSFINCS ]; then\
		echo "----SFINCS OPTION----";\
		echo $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.vboot_tokamak_sfincs;\
		cd $(STELLOPT_TEST_DIR)/vboot_tokamak_sfincs; $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.vboot_tokamak_sfincs;\
		echo "";\
	fi
	@echo "Tesing of STELLOPT complete"
	@echo ""

help:
	@echo "TESTING OF STELLOPT BUILDS"
	@echo "   ENVIRONMENT VARIABLES (make.inc)"
	@echo "      MPI_RUN: " $(MPI_RUN)
	@echo "      MPI_RUN_OPTS: " $(MPI_RUN_OPTS)
	@echo ""
	@echo "   COMMANDS"
	@echo "      test_vmec     Tesing of VMEC"
	@echo "      test_diagno   Tesing of DIAGNO"
	@echo "      test_stellopt Tesing of STELLOPT"
