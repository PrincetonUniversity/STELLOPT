#Type "make release(debug)" at command line to make release(debug) version  
#DO NOT alter the directory structure or move this makefile (must be in first subdirectory below top_level_dir)
  MAKECMDGOALS = xstelloptv2
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
	@echo ""
	@echo $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xvmec2000 input.DIII-D
	@cd $(VMEC_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xvmec2000 input.DIII-D
	@cd $(VMEC_TEST_DIR); ./compare_DIII-D.py
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
	$(MAKE) test_stellopt_basic
	$(MAKE) test_stellopt_iota_lmdif
	$(MAKE) test_stellopt_tok_r0_delta
	$(MAKE) test_stellopt_tok_r0_rho
	$(MAKE) test_stellopt_qhs_lmdif
	#$(MAKE) test_stellopt_beams3d
	$(MAKE) test_stellopt_vboot_tokamak_bootsj
	@if [ $LCOILOPT ]; then\
		@echo "----COILOPT OPTION----"
		$(MAKE) test_stellopt_coilopt
	fi
	@if [ $LGENE ]; then\
		@echo "----GENE OPTION----"
		$(MAKE) test_stellopt_gene
	fi
	#@if [ $LREGCOIL ]; then\
	#	@echo "----REGCOIL OPTION----"
	#	$(MAKE) test_stellopt_regcoil
	#fi
	@if [ $LSFINCS ]; then\
		@echo "----SFINCS OPTION----"
		$(MAKE) test_stellopt_sfincs
	fi
	@if [ $LTERPSICHORE ]; then\
		@echo "----TERPSICHORE OPTION----"
		$(MAKE) test_stellopt_terpsichore
	fi
	#@if [ $LTRAVIS ]; then\
	#	@echo "----TRAVIS OPTION----"
	#	$(MAKE) test_stellopt_travis
	#fi
	@echo "Tesing of STELLOPT complete"
	@echo ""

test_stellopt_basic:
	@echo $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.BASIC
	@cd $(STELLOPT_TEST_DIR)/BASIC; $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.BASIC; rm -f *_opt*
	@echo ""

test_stellopt_beams3d:
	@echo $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.BEAMS3D
	@cd $(STELLOPT_TEST_DIR)/BEAMS3D; $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.BEAMS3D; rm -f *_opt*
	@echo ""

test_stellopt_iota_lmdif:
	@echo $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.IOTA_LMDIF
	@cd $(STELLOPT_TEST_DIR)/IOTA_LMDIF; $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.IOTA_LMDIF; rm -f *_opt*
	@echo ""

test_stellopt_tok_r0_delta:
	@echo $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.TOK_R0_DELTA
	@cd $(STELLOPT_TEST_DIR)/TOK_R0_DELTA; $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.TOK_R0_DELTA; rm -f *_opt*
	@echo ""

test_stellpt_tok_r0_rho:
	@echo $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.TOK_R0_RHO
	@cd $(STELLOPT_TEST_DIR)/TOK_R0_RHO; $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.TOK_R0_RHO; rm -f *_opt*
	@echo ""

test_stellopt_qhs_lmdif:
	@echo $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.QHS_LMDIF
	@cd $(STELLOPT_TEST_DIR)/QHS_LMDIF; $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.QHS_LMDIF; rm -f *_opt*
	@echo ""

test_stellopt_vboot_tokamak_bootsj:
	@echo $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.vboot_tokamak_bootsj
	@cd $(STELLOPT_TEST_DIR)/vboot_tokamak_bootsj; $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.vboot_tokamak_bootsj; rm -f *_opt*
	@echo ""

test_stellopt_coilopt:
	@echo $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.COILOPTPP
	@cd $(STELLOPT_TEST_DIR)/COILOPTPP; $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.COILOPTPP; rm -f *_opt*
	@echo ""

test_stellopt_gene:
	@echo $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.GENE;
	@cd $(STELLOPT_TEST_DIR)/GENE; $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.GENE; rm -f *_opt*
	@echo ""

test_stellopt_sfincs:
	@echo $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.vboot_tokamak_sfincs;
	@cd $(STELLOPT_TEST_DIR)/vboot_tokamak_sfincs; $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.vboot_tokamak_sfincs; rm -f *_opt*
	@echo ""

test_stellopt_terpsichore:
	@echo $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.TERPSICHORE
	@cd $(STELLOPT_TEST_DIR)/TERPSICHORE; $(MPI_RUN) $(MPI_RUN_OPTS) $(MYHOME)/xstelloptv2 input.TERPSICHORE; rm -f *_opt*
	@echo ""

clean_test:
	@cd $(VMEC_TEST_DIR); rm -f dcon* jxbout* mercier* threed1* parvmecinfo* wout*
	@cd $(DIAGNO_TEST_DIR); rm -f diagno_bench* diagno_flux* diagno_surf* test_loops_j* test_loops_b*
	@$(MAKE) SUBDIR=$(STELLOPT_TEST_DIR)/BASIC stellopt_cleanup_dir
	@$(MAKE) SUBDIR=$(STELLOPT_TEST_DIR)/BEAMS3D stellopt_cleanup_dir
	@$(MAKE) SUBDIR=$(STELLOPT_TEST_DIR)/COILOPTPP stellopt_cleanup_dir
	@$(MAKE) SUBDIR=$(STELLOPT_TEST_DIR)/GENE stellopt_cleanup_dir
	@$(MAKE) SUBDIR=$(STELLOPT_TEST_DIR)/IOTA_LMDIF stellopt_cleanup_dir
	@$(MAKE) SUBDIR=$(STELLOPT_TEST_DIR)/ONE_ITER stellopt_cleanup_dir
	@$(MAKE) SUBDIR=$(STELLOPT_TEST_DIR)/QHS_LMDIF stellopt_cleanup_dir
	@$(MAKE) SUBDIR=$(STELLOPT_TEST_DIR)/TERPSICHORE stellopt_cleanup_dir
	@$(MAKE) SUBDIR=$(STELLOPT_TEST_DIR)/BASIC stellopt_cleanup_dir
	@$(MAKE) SUBDIR=$(STELLOPT_TEST_DIR)/TOK_R0_DELTA stellopt_cleanup_dir
	@$(MAKE) SUBDIR=$(STELLOPT_TEST_DIR)/TOK_R0_RHO stellopt_cleanup_dir
	@$(MAKE) SUBDIR=$(STELLOPT_TEST_DIR)/vboot_tokamak_sfincs stellopt_cleanup_dir

stellopt_cleanup_dir:
	@cd $(SUBDIR); rm -f answers* auto* beams3d* b_norm_eq_* b_norm_init_* b_norm_final_* bnorm.* boot_fit.* boozmn* checkpoint_* \
		coilopt_params.* coils.* coil_spline*.out curv* dcon* deltamn.* diagno_flux.* diagno_bth.* diagno_seg.* \
		dkes* dprof* fevals* field* fort.* gade* gist* eigenvalues* input.*.***** input_dkes.* jacobian.* jprof.* \
		jxb* jBbs* LMsave* log_gene.* map.* map_plane.* mgrid_*.nc mercier* mom* neo* nrg* opt* parameters_* \
		parvmecinfo* results* rhomn* RZ* stellopt.* tprof* wout* temp_input.* terpsichore_eq.* terpsichore_16.* \
		terpsichore_17.* terpsichore_19.* terpsichore_22.* terpsichore_23.* threed1.* tv_out_* txport* xvec.* uv_* \
		var_* vsp*
	

help:
	@echo "TESTING OF STELLOPT BUILDS"
	@echo "   ENVIRONMENT VARIABLES (make.inc)"
	@echo "      MPI_RUN:      " $(MPI_RUN)
	@echo "      MPI_RUN_OPTS: " $(MPI_RUN_OPTS)
	@echo ""
	@echo "      LCOILOPT:     " $(LCOILOPT)
	@echo "      LGENE:        " $(LGENE)
	@echo "      LREGCOIL:     " $(LREGCOIL)
	@echo "      LSINFCS:      " $(LSFINCS)
	@echo "      LTERPSICHORE: " $(LTERPSICHORE)
	@echo "      LTRAVIS:      " $(LTRAVIS)
	@echo ""
	@echo "   COMMANDS"
	@echo "      test_vmec     Tesing of VMEC"
	@echo "      test_diagno   Tesing of DIAGNO"
	@echo "      test_stellopt Tesing of STELLOPT"
	@echo "      clean_test    Remove files generated by tests"
