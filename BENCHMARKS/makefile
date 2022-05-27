#Basic files
  include ../make.inc

# BINARY DIRECTORIES
  LIB_DIR = $(MYHOME)/libstell_dir

# Setup Subdirectories
  BEAMS3D_TEST_DIR = BEAMS3D_TEST
  DIAGNO_TEST_DIR = DIAGNO_TEST
  FIELDLINES_TEST_DIR = FIELDLINES_TEST
  STELLOPT_TEST_DIR = STELLOPT_TEST
  VMEC_TEST_DIR = VMEC_TEST

# We do the following to handle Testing on Large Systems
  MPI_RUN_OPTS_SM ?= $(MPI_RUN_OPTS)
  MPI_RUN_OPTS_MD ?= $(MPI_RUN_OPTS)
  MPI_RUN_OPTS_LG ?= $(MPI_RUN_OPTS)

.PHONY: test_vmec

test_vmec: 
	@echo "Beginning testing of VMEC2000"
	@cd $(VMEC_TEST_DIR); rm -f jxbout* mercier* threed1* wout*; 
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xvmec2000 input.DSHAPE
	@cd $(VMEC_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xvmec2000 input.DSHAPE
	@cd $(VMEC_TEST_DIR) && ./compare_DSHAPE.py
	@echo ""
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xvmec2000 input.SOLOVEV
	@cd $(VMEC_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xvmec2000 input.SOLOVEV
	@cd $(VMEC_TEST_DIR); ./compare_SOLOVEV.py
	@echo ""
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xvmec2000 input.HELIOTRON
	@cd $(VMEC_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xvmec2000 input.HELIOTRON
	@cd $(VMEC_TEST_DIR); ./compare_HELIOTRON.py
	@echo ""
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xvmec2000 input.ATF
	@cd $(VMEC_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xvmec2000 input.ATF
	@cd $(VMEC_TEST_DIR); ./compare_ATF.py
	@echo ""
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xvmec2000 input.QAS
	@cd $(VMEC_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xvmec2000 input.QAS
	@cd $(VMEC_TEST_DIR); ./compare_QAS.py
	@echo ""
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xvmec2000 input.DIII-D
	@cd $(VMEC_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xvmec2000 input.DIII-D
	@cd $(VMEC_TEST_DIR); ./compare_DIII-D.py
	@echo "Tesing of VMEC2000 complete"

test_diagno:
	@echo "Beginning testing of DIAGNO"
	@cd $(DIAGNO_TEST_DIR); rm -f diagno_bench.* diagno_flux.* diagno_surf.* test_loops_b.* test_loops_j.*; 
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_LG) $(MYHOME)/xdiagno -vmec bigtok -bench
	@cd $(DIAGNO_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS_LG) $(MYHOME)/xdiagno -bench -vmec bigtok
	@cd $(DIAGNO_TEST_DIR); ./compare_bigtok.py
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_LG) $(MYHOME)/xdiagno -vmec bigtok_phi1 -bench
	@cd $(DIAGNO_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS_LG) $(MYHOME)/xdiagno -bench -vmec bigtok_phi1
	@cd $(DIAGNO_TEST_DIR); ./compare_bigtok_phi1.py
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_LG) $(MYHOME)/xdiagno -vmec a3tok -bench
	@cd $(DIAGNO_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS_LG) $(MYHOME)/xdiagno -bench -vmec a3tok
	@cd $(DIAGNO_TEST_DIR); ./compare_a3tok.py
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_LG) $(MYHOME)/xdiagno -vmec DIIID_m24n0s99_nfp1 -bench
	@cd $(DIAGNO_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS_LG) $(MYHOME)/xdiagno -bench -vmec DIIID_m24n0s99_nfp1
	@cd $(DIAGNO_TEST_DIR); ./compare_DIIID_m24n0s99_nfp1.py
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_LG) $(MYHOME)/xdiagno -vmec DIIID_m24n0s99_nfp10 -bench
	@cd $(DIAGNO_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS_LG) $(MYHOME)/xdiagno -bench -vmec DIIID_m24n0s99_nfp10
	@cd $(DIAGNO_TEST_DIR); ./compare_DIIID_m24n0s99_nfp10.py
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_LG) $(MYHOME)/xdiagno -vmec DIIID_m20n0s128_nfp1_lasym -bench
	@cd $(DIAGNO_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS_LG) $(MYHOME)/xdiagno -bench -vmec DIIID_m20n0s128_nfp1_lasym
	@cd $(DIAGNO_TEST_DIR); ./compare_DIIID_m20n0s128_nfp1_lasym.py
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_LG) $(MYHOME)/xdiagno -vmec iter -bench
	@cd $(DIAGNO_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS_LG) $(MYHOME)/xdiagno -bench -vmec iter
	@cd $(DIAGNO_TEST_DIR); ./compare_iter.py
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_LG) $(MYHOME)/xdiagno -vmec lhd -bench
	@cd $(DIAGNO_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS_LG) $(MYHOME)/xdiagno -bench -vmec lhd
	@cd $(DIAGNO_TEST_DIR); ./compare_lhd.py
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_LG) $(MYHOME)/xdiagno -vmec hsx -bench
	@cd $(DIAGNO_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS_LG) $(MYHOME)/xdiagno -bench -vmec hsx
	@cd $(DIAGNO_TEST_DIR); ./compare_hsx.py
	@echo ""

test_fieldlines:
	@echo "Beginning testing of FIELDLINES"
	$(MAKE) test_fieldlines_full
	$(MAKE) test_fieldlines_coll
	@echo "Tesing of FIELDLINES complete"
	@echo ""

test_beams3d:
	@echo "Beginning testing of BEAMS3D"
	$(MAKE) test_beams3d_basic
	$(MAKE) test_beams3d_er
	$(MAKE) test_beams3d_loss
	$(MAKE) test_beams3d_slow
	$(MAKE) test_beams3d_restart
	$(MAKE) test_beams3d_eqdsk
	$(MAKE) test_beams3d_fusion
	$(MAKE) test_beams3d_multiion
	@echo "Tesing of BEAMS3D complete"
	@echo ""

test_beams3d_basic:
	@cd $(BEAMS3D_TEST_DIR); rm -f *.h5 *.txt; 
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xbeams3d -vmec ORBITS -plasma
	@cd $(BEAMS3D_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xbeams3d -vmec ORBITS -plasma
	@cd $(BEAMS3D_TEST_DIR); ./test_ORBITS.py
	@echo ""

test_beams3d_er:
	@cd $(BEAMS3D_TEST_DIR); rm -f *.h5 *.txt; 
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xbeams3d -vmec ORBITS_er -plasma
	@cd $(BEAMS3D_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xbeams3d -vmec ORBITS_er -plasma
	@cd $(BEAMS3D_TEST_DIR); ./test_ORBITS_er.py
	@echo ""
	@echo ""

test_beams3d_loss:
	@cd $(BEAMS3D_TEST_DIR); rm -f *.h5 *.txt; 
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xbeams3d -vmec ORBITS_loss -plasma
	@cd $(BEAMS3D_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xbeams3d -vmec ORBITS_loss -plasma
	@cd $(BEAMS3D_TEST_DIR); ./test_ORBITS_loss.py
	@echo ""

test_beams3d_slow:
	@cd $(BEAMS3D_TEST_DIR); rm -f *.h5 *.txt; 
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xbeams3d -vmec ORBITS_slow -plasma -collisions
	@cd $(BEAMS3D_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xbeams3d -vmec ORBITS_slow -plasma -collisions
	@cd $(BEAMS3D_TEST_DIR); ./test_ORBITS_slow.py
	@echo ""

test_beams3d_depo:
	@cd $(BEAMS3D_TEST_DIR); rm -f *.h5 *.txt; 
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xbeams3d -vmec ORBITS_depo -plasma -depo -suzuki
	@cd $(BEAMS3D_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xbeams3d -vmec ORBITS_depo -plasma -depo -suzuki
	@cd $(BEAMS3D_TEST_DIR); ./test_ORBITS_depo.py
	@echo ""

test_beams3d_restart:
	$(MAKE) test_beams3d_depo
	@cd $(BEAMS3D_TEST_DIR);
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xbeams3d -vmec ORBITS_restart -plasma -restart beams3d_ORBITS_depo.h5
	@cd $(BEAMS3D_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xbeams3d -vmec ORBITS_restart -plasma -restart beams3d_ORBITS_depo.h5
	@cd $(BEAMS3D_TEST_DIR); ./test_ORBITS_restart.py
	@echo ""

test_beams3d_dist:
	@cd $(BEAMS3D_TEST_DIR); rm -f *.h5 *.txt; 
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xbeams3d -vmec ORBITS_dist -plasma -collisions
	@cd $(BEAMS3D_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xbeams3d -vmec ORBITS_dist -plasma -collisions
	@cd $(BEAMS3D_TEST_DIR); ./test_ORBITS_dist.py
	@echo ""

test_beams3d_fusion:
	@cd $(BEAMS3D_TEST_DIR); rm -f *.h5 *.txt; 
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xbeams3d -vmec ORBITS_fusion -plasma -fusion -collisions -depo
	@cd $(BEAMS3D_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xbeams3d -vmec ORBITS_fusion -plasma -fusion -collisions -depo
	@cd $(BEAMS3D_TEST_DIR); ./test_ORBITS_fusion.py
	@echo ""

test_beams3d_eqdsk:
	@cd $(BEAMS3D_TEST_DIR); rm -f *.h5 *.txt; 
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xbeams3d -eqdsk ORBITS_eqdsk g164723.03059 -depo -suzuki
	@cd $(BEAMS3D_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xbeams3d -eqdsk ORBITS_eqdsk g164723.03059 -depo -suzuki
	@cd $(BEAMS3D_TEST_DIR); ./test_ORBITS_eqdsk.py
	@echo ""

test_beams3d_multiion:
	@cd $(BEAMS3D_TEST_DIR); rm -f *.h5 *.txt; 
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xbeams3d -vmec ORBITS_multiion -plasma -depo -suzuki
	@cd $(BEAMS3D_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xbeams3d -vmec ORBITS_multiion -plasma -depo -suzuki
	@cd $(BEAMS3D_TEST_DIR); ./test_ORBITS_multiion.py
	@echo ""

test_beams3d_fullorbit:
	@cd $(BEAMS3D_TEST_DIR); rm -f *.h5 *.txt; 
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xbeams3d -vmec ORBITS_fullorbit -plasma
	@cd $(BEAMS3D_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xbeams3d -vmec ORBITS_fullorbit -plasma
	#@cd $(BEAMS3D_TEST_DIR); ./test_beams3d_fullorbit.py
	@echo ""

test_stellopt:
	@echo "Beginning testing of STELLOPT"
	$(MAKE) test_stellopt_lmdif
	$(MAKE) test_stellopt_gade
	$(MAKE) test_stellopt_basic
	$(MAKE) test_stellopt_iota_lmdif
	$(MAKE) test_stellopt_tok_r0_delta
	$(MAKE) test_stellopt_tok_r0_rho
	$(MAKE) test_stellopt_qhs_lmdif
	$(MAKE) test_stellopt_beams3d
	$(MAKE) test_stellopt_recon_tok
	$(MAKE) test_stellopt_vboot_tokamak_bootsj
	@if [ $LCOILOPT ]; then\
		echo "----COILOPT OPTION----"; \
		$(MAKE) test_stellopt_coilopt;\
	fi
	@if [ $LGENE ]; then\
		@echo "----GENE OPTION----";\
		$(MAKE) test_stellopt_gene;\
	fi
	@if [ $LREGCOIL ]; then\
		@echo "----REGCOIL OPTION----";\
		$(MAKE) test_stellopt_regcoil;\
	fi
	@if [ $LSFINCS ]; then\
		@echo "----SFINCS OPTION----";\
		$(MAKE) test_stellopt_sfincs;\
	fi
	@if [ $LTERPSICHORE ]; then\
		@echo "----TERPSICHORE OPTION----";\
		$(MAKE) test_stellopt_terpsichore;\
	fi
	@if [ $LTRAVIS ]; then\
		@echo "----TRAVIS OPTION----";\
		$(MAKE) test_stellopt_travis;\
	fi
	@echo "Tesing of STELLOPT complete"
	@echo ""

test_stellopt_lmdif:
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_MD) $(MYHOME)/xstelloptv2 input.LMDIF_TEST
	@cd $(STELLOPT_TEST_DIR)/LMDIF_TEST; $(MPI_RUN) $(MPI_RUN_OPTS_MD) $(MYHOME)/xstelloptv2 input.LMDIF_TEST; rm -f *_opt*
	@cd $(STELLOPT_TEST_DIR)/LMDIF_TEST && ./compare_LMDIF_TEST.py
	@echo ""

test_stellopt_gade:
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_MD) $(MYHOME)/xstelloptv2 input.GADE_TEST
	@cd $(STELLOPT_TEST_DIR)/GADE_TEST; $(MPI_RUN) $(MPI_RUN_OPTS_MD) $(MYHOME)/xstelloptv2 input.GADE_TEST; rm -f *_opt*
	@cd $(STELLOPT_TEST_DIR)/GADE_TEST && ./compare_GADE_TEST.py
	@echo ""

test_stellopt_basic:
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_MD) $(MYHOME)/xstelloptv2 input.BASIC
	@cd $(STELLOPT_TEST_DIR)/BASIC; $(MPI_RUN) $(MPI_RUN_OPTS_MD) $(MYHOME)/xstelloptv2 input.BASIC; rm -f *_opt*
	@cd $(STELLOPT_TEST_DIR)/BASIC && ./compare_BASIC.py
	@echo ""

test_stellopt_beams3d:
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_LG) $(MYHOME)/xstelloptv2 input.BEAMS3D
	@cd $(STELLOPT_TEST_DIR)/BEAMS3D; $(MPI_RUN) $(MPI_RUN_OPTS_LG) $(MYHOME)/xstelloptv2 input.BEAMS3D; rm -f *_opt*
	@cd $(STELLOPT_TEST_DIR)/BEAMS3D && ./compare_STELLOPT_BEAMS3D.py
	@echo ""

test_stellopt_iota_lmdif:
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_MD) $(MYHOME)/xstelloptv2 input.IOTA_LMDIF
	@cd $(STELLOPT_TEST_DIR)/IOTA_LMDIF; $(MPI_RUN) $(MPI_RUN_OPTS_MD) $(MYHOME)/xstelloptv2 input.IOTA_LMDIF; rm -f *_opt*
	@cd $(STELLOPT_TEST_DIR)/IOTA_LMDIF && ./compare_IOTA_LMDIF.py
	@echo ""

test_stellopt_recon_tok:
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_MD) $(MYHOME)/xstelloptv2 input.RECON_TOK
	@cd $(STELLOPT_TEST_DIR)/RECON_TOK; $(MPI_RUN) $(MPI_RUN_OPTS_MD) $(MYHOME)/xstelloptv2 input.RECON_TOK; rm -f *_opt*
	@cd $(STELLOPT_TEST_DIR)/RECON_TOK && ./compare_RECON_TOK.py
	@echo ""

test_stellopt_tok_r0_delta:
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_MD) $(MYHOME)/xstelloptv2 input.TOK_R0_DELTA
	@cd $(STELLOPT_TEST_DIR)/TOK_R0_DELTA; $(MPI_RUN) $(MPI_RUN_OPTS_MD) $(MYHOME)/xstelloptv2 input.TOK_R0_DELTA; rm -f *_opt*
	@cd $(STELLOPT_TEST_DIR)/TOK_R0_DELTA && ./compare_TOK_R0_DELTA.py
	@echo ""

test_stellopt_tok_r0_rho:
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_MD) $(MYHOME)/xstelloptv2 input.TOK_R0_RHO
	@cd $(STELLOPT_TEST_DIR)/TOK_R0_RHO; $(MPI_RUN) $(MPI_RUN_OPTS_MD) $(MYHOME)/xstelloptv2 input.TOK_R0_RHO; rm -f *_opt*
	@cd $(STELLOPT_TEST_DIR)/TOK_R0_RHO && ./compare_TOK_R0_RHO.py
	@echo ""

test_stellopt_qhs_lmdif:
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_MD) $(MYHOME)/xstelloptv2 input.QHS_LMDIF
	@cd $(STELLOPT_TEST_DIR)/QHS_LMDIF; $(MPI_RUN) $(MPI_RUN_OPTS_MD) $(MYHOME)/xstelloptv2 input.QHS_LMDIF; rm -f *_opt*
	@cd $(STELLOPT_TEST_DIR)/QHS_LMDIF && ./compare_QHS_LMDIF.py
	@echo ""

test_stellopt_vboot_tokamak_bootsj:
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_MD) $(MYHOME)/xstelloptv2 input.vboot_tokamak_bootsj
	@cd $(STELLOPT_TEST_DIR)/vboot_tokamak_bootsj; $(MPI_RUN) $(MPI_RUN_OPTS_MD) $(MYHOME)/xstelloptv2 input.vboot_tokamak_bootsj; rm -f *_opt*
	@echo ""

test_stellopt_coilopt:
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_LG) $(MYHOME)/xstelloptv2 input.COILOPTPP
	@cd $(STELLOPT_TEST_DIR)/COILOPTPP; $(MPI_RUN) $(MPI_RUN_OPTS_LG) $(MYHOME)/xstelloptv2 input.COILOPTPP; rm -f *_opt*
	@echo ""

test_stellopt_gene:
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_LG) $(MYHOME)/xstelloptv2 input.GENE;
	@cd $(STELLOPT_TEST_DIR)/GENE; $(MPI_RUN) $(MPI_RUN_OPTS_LG) $(MYHOME)/xstelloptv2 input.GENE; rm -f *_opt*
	@echo ""

test_stellopt_sfincs:
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_MD) $(MYHOME)/xstelloptv2 input.vboot_tokamak_sfincs;
	@cd $(STELLOPT_TEST_DIR)/vboot_tokamak_sfincs; $(MPI_RUN) $(MPI_RUN_OPTS_MD) $(MYHOME)/xstelloptv2 input.vboot_tokamak_sfincs; rm -f *_opt*
	@echo ""

test_stellopt_terpsichore:
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_MD) $(MYHOME)/xstelloptv2 input.TERPSICHORE
	@cd $(STELLOPT_TEST_DIR)/TERPSICHORE; $(MPI_RUN) $(MPI_RUN_OPTS_MD) $(MYHOME)/xstelloptv2 input.TERPSICHORE; rm -f *_opt*
	@echo ""

test_stellopt_travis:
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xstelloptv2 input.TRAVIS
	@cd $(STELLOPT_TEST_DIR)/TRAVIS; $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xstelloptv2 input.TRAVIS; rm -f *_opt*
	@cd $(STELLOPT_TEST_DIR)/TRAVIS && ./compare_TRAVIS.py
	@echo ""

test_fieldlines_full:
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xfieldlines -vmec NCSX_s1 -coil coils.NCSX -vac -full
	@cd $(FIELDLINES_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS_SM) $(MYHOME)/xfieldlines -vmec NCSX_s1 -vac -coil coils.NCSX -vac -full
	@cd $(FIELDLINES_TEST_DIR); ./compare_NCSX_s1.py
	@echo ""

test_fieldlines_coll: FIELDLINES_TEST/fieldlines_NCSX_s1.h5
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_MD) $(MYHOME)/xfieldlines -vmec NCSX_s1_coll -coil coils.NCSX -vac -vessel NCSX_wall_trimesh.dat -hitonly -field_start fieldlines_NCSX_s1.h5 114
	@cd $(FIELDLINES_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS_MD) $(MYHOME)/xfieldlines -vmec NCSX_s1_coll -vac -coil coils.NCSX -vac -vessel NCSX_wall_trimesh.dat -hitonly -field_start fieldlines_NCSX_s1.h5 114
	@echo ""

test_fieldlines_backflow: FIELDLINES_TEST/fieldlines_NCSX_s1.h5
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_MD) $(MYHOME)/xfieldlines -vmec NCSX_s1_coll -coil coils.NCSX -vac -vessel NCSX_wall_trimesh.dat -hitonly -field_start fieldlines_NCSX_s1.h5 114 -backflow
	@cd $(FIELDLINES_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS_MD) $(MYHOME)/xfieldlines -vmec NCSX_s1_coll -vac -coil coils.NCSX -vac -vessel NCSX_wall_trimesh.dat -hitonly -field_start fieldlines_NCSX_s1.h5 114 -backflow
	@echo ""

test_fieldlines_coll_acc: FIELDLINES_TEST/fieldlines_NCSX_s1.h5
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_MD) $(MYHOME)/xfieldlines -vmec NCSX_s1_coll -coil coils.NCSX -vac -vessel NCSX_wall_trimesh_acc0.25v4.dat -hitonly -field_start fieldlines_NCSX_s1.h5 114
	@cd $(FIELDLINES_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS_MD) $(MYHOME)/xfieldlines -vmec NCSX_s1_coll -vac -coil coils.NCSX -vac -vessel NCSX_wall_trimesh_acc0.25v4.dat -hitonly -field_start fieldlines_NCSX_s1.h5 114
	@echo ""

test_fieldlines_coll_ef: FIELDLINES_TEST/fieldlines_NCSX_s1.h5
	@echo $(MPI_RUN) $(MPI_RUN_OPTS_MD) $(MYHOME)/xfieldlines -vmec NCSX_s1_coll_ef -coil coils.NCSX_nfp1 -vac -vessel NCSX_wall_trimesh.dat -hitonly -field_start fieldlines_NCSX_s1.h5 114
	@cd $(FIELDLINES_TEST_DIR); $(MPI_RUN) $(MPI_RUN_OPTS_MD) $(MYHOME)/xfieldlines -vmec NCSX_s1_coll_ef -vac -coil coils.NCSX_nfp1 -vac -vessel NCSX_wall_trimesh.dat -hitonly -field_start fieldlines_NCSX_s1.h5 114
	@echo ""

FIELDLINES_TEST/fieldlines_NCSX_s1.h5: 
	@$(MAKE) test_fieldlines_full

clean_test:
	@cd $(VMEC_TEST_DIR); rm -f dcon* jxbout* mercier* threed1* parvmecinfo* wout*
	@cd $(DIAGNO_TEST_DIR); rm -f diagno_bench* diagno_flux* diagno_surf* test_loops_j* test_loops_b*
	@cd $(FIELDLINES_TEST_DIR); rm -f fieldlines_*.h5
	@$(MAKE) SUBDIR=$(STELLOPT_TEST_DIR)/LMDIF_TEST stellopt_cleanup_dir
	@$(MAKE) SUBDIR=$(STELLOPT_TEST_DIR)/GADE_TEST stellopt_cleanup_dir
	@$(MAKE) SUBDIR=$(STELLOPT_TEST_DIR)/BASIC stellopt_cleanup_dir
	@$(MAKE) SUBDIR=$(STELLOPT_TEST_DIR)/BEAMS3D stellopt_cleanup_dir
	@$(MAKE) SUBDIR=$(STELLOPT_TEST_DIR)/COILOPTPP stellopt_cleanup_dir
	@$(MAKE) SUBDIR=$(STELLOPT_TEST_DIR)/GENE stellopt_cleanup_dir
	@$(MAKE) SUBDIR=$(STELLOPT_TEST_DIR)/IOTA_LMDIF stellopt_cleanup_dir
	@$(MAKE) SUBDIR=$(STELLOPT_TEST_DIR)/ONE_ITER stellopt_cleanup_dir
	@$(MAKE) SUBDIR=$(STELLOPT_TEST_DIR)/QHS_LMDIF stellopt_cleanup_dir
	@$(MAKE) SUBDIR=$(STELLOPT_TEST_DIR)/RECON_TOK stellopt_cleanup_dir
	@$(MAKE) SUBDIR=$(STELLOPT_TEST_DIR)/TERPSICHORE stellopt_cleanup_dir
	@$(MAKE) SUBDIR=$(STELLOPT_TEST_DIR)/TOK_R0_DELTA stellopt_cleanup_dir
	@$(MAKE) SUBDIR=$(STELLOPT_TEST_DIR)/TOK_R0_RHO stellopt_cleanup_dir
	@$(MAKE) SUBDIR=$(STELLOPT_TEST_DIR)/TRAVIS stellopt_cleanup_dir
	@$(MAKE) SUBDIR=$(STELLOPT_TEST_DIR)/vboot_tokamak_sfincs stellopt_cleanup_dir

stellopt_cleanup_dir:
	@cd $(SUBDIR); rm -f answers* auto* beams3d* b_norm_eq_* b_norm_init_* b_norm_final_* bnorm.* boot_fit.* boozmn* checkpoint_* \
		coilopt_params.* coils.* coil_spline*.out curv* dcon* deltamn.* diagno_flux.* diagno_bth.* diagno_seg.* \
		dkes* dprof* fevals* field* fort.* gade* gist* eigenvalues* input.*_min input.*.***** input_dkes.* jacobian.* jprof.* \
		jxb* jBbs* LMsave* log_gene.* map.* map_plane.* mgrid_*.nc mercier* mom* neo* nrg* opt* parameters_* \
		parvmecinfo* results* rhomn* RZ* stellopt.* tprof* wout* temp_input.* terpsichore_eq.* terpsichore_16.* \
		terpsichore_17.* terpsichore_19.* terpsichore_22.* terpsichore_23.* threed1.* tv_out_* txport* xvec.* uv_* \
		var_* vsp*
	

help:
	@echo "TESTING OF STELLOPT BUILDS"
	@echo "   ENVIRONMENT VARIABLES (make.inc)"
	@echo "      MPI_RUN:         " $(MPI_RUN)
	@echo "      MPI_RUN_OPTS:    " $(MPI_RUN_OPTS)
	@echo "      MPI_RUN_OPTS_SM: " $(MPI_RUN_OPTS_SM)
	@echo "      MPI_RUN_OPTS_MD: " $(MPI_RUN_OPTS_MD)
	@echo "      MPI_RUN_OPTS_LG: " $(MPI_RUN_OPTS_LG)
	@echo ""
	@echo "      LCOILOPT:        " $(LCOILOPT)
	@echo "      LGENE:           " $(LGENE)
	@echo "      LREGCOIL:        " $(LREGCOIL)
	@echo "      LSINFCS:         " $(LSFINCS)
	@echo "      LTERPSICHORE:    " $(LTERPSICHORE)
	@echo "      LTRAVIS:         " $(LTRAVIS)
	@echo ""
	@echo "   COMMANDS"
	@echo "      test_beams3d    Tesing of BEAMS3D"
	@echo "      test_diagno     Tesing of DIAGNO"
	@echo "      test_fieldlines Tesing of FIELDLINES"
	@echo "      test_stellopt   Tesing of STELLOPT"
	@echo "      test_vmec       Tesing of VMEC"
	@echo "      clean_test      Remove files generated by tests"
