###############################################################
# Variables

PFUNIT_TEST_DIR = $(GK_HEAD_DIR)/tests/unit
PFUNIT_NPROC = 4

# Default major version of pFUnit
# Versions 3 and 4 require different flags, etc.
PFUNIT_VERSION_MAJOR ?= 4
# Installation directory of pFUnit
# Default is the directory where `scripts/build_pfunit` puts it
PFUNIT_DIR ?= $(GK_HEAD_DIR)/externals/pFUnit_build/install/PFUNIT-4.2

PFUNIT_LOG:=pfunit.log
PFUNIT_ERROR:=pfunit.error

CPPFLAGS += -DPFUNIT_VERSION_MAJOR=$(PFUNIT_VERSION_MAJOR)

ifeq ($(PFUNIT_VERSION_MAJOR),4)
  sinclude $(PFUNIT_DIR)/include/PFUNIT.mk
  PFUNIT_MOD := pfunit
  PFUNIT_PREPROCESSOR ?= $(PFUNIT_BIN_DIR)/funitproc
  PFUNIT_INC := -I$(PFUNIT_INCLUDE_DIR)
  PFUNIT_LIB := -L$(PFUNIT_TOP_DIR)/lib -lpfunit -lfunit $(FARGPARSE_LIBRARIES)
  F90FLAGS += $(PFUNIT_EXTRA_FFLAGS)
else
  PFUNIT_MOD := pfunit_mod
  PFUNIT_PREPROCESSOR ?= $(PFUNIT_DIR)/bin/pFUnitParser.py
  PFUNIT_INC := -I$(PFUNIT_DIR)/mod
  PFUNIT_LIB := -L$(PFUNIT_DIR)/lib -lpfunit
endif

CPPFLAGS += -DPFUNIT_MOD=$(PFUNIT_MOD)

PFUNIT_TEST_PATH :=.

STELLA_OBJS := $(filter-out stella.o, $(stella_mod))

PFUNIT_TESTS:= $(shell find $(PFUNIT_TEST_DIR) -name '[^.\#]*.pf')

PFUNIT_OBJS := $(PFUNIT_TESTS:.pf=.o)
PFUNIT_INTERMEDIATES := $(PFUNIT_TESTS:.pf=.F90)
PFUNIT_DRIVER := $(PFUNIT_DIR)/include/driver.F90
PFUNIT_TEST_LIST := $(PFUNIT_TEST_DIR)/testSuites.inc

.INTERMEDIATE: $(PFUNIT_INTERMEDIATES)

###############################################################
# pFUnit library

.PHONY: build-pfunit-library check-pfunit-install clean-pfunit-library

build-pfunit-library:
	@echo "Building pFUnit in default location: $(GK_HEAD_DIR)/externals/pFUnit_build/install"
ifdef USE_MPI
	@scripts/build_pfunit -v ${PFUNIT_VERSION_MAJOR}
else
	@scripts/build_pfunit -v ${PFUNIT_VERSION_MAJOR} -m
endif

check-pfunit-install:
ifeq ("$(wildcard $(PFUNIT_PREPROCESSOR))","")
	$(error Could not find $(PFUNIT_PREPROCESSOR). Make sure pFUnit is installed and PFUNIT_DIR ($(PFUNIT_DIR)) is set to the current location)
endif

clean-pfunit-library:
	@echo "Removing pFUnit library"
	@scripts/build_pfunit -r

distclean: clean-pfunit-library

###############################################################
# Test suite

%.o : %.pf
	$(PFUNIT_PREPROCESSOR) $< $(basename $<).F90
ifdef ONE_STEP_PP
	$(QUIETSYM)$(FC) $(FPPFLAGS) $(F90FLAGS) -c $(basename $<).F90 $(PFUNIT_INC) -I$(GK_HEAD_DIR) -o $@
else
	$(QUIETSYM)$(CPP) $(FPPFLAGS) $(basename $<).F90 -o $(basename $<).f90
	$(QUIETSYM)$(FC) $(F90FLAGS) -c $(basename $<).f90 $(PFUNIT_INC) -I$(GK_HEAD_DIR) -o $@
endif

$(PFUNIT_DRIVER): $(PFUNIT_TEST_LIST)
	@touch $(PFUNIT_DRIVER)

$(PFUNIT_OBJS): $(STELLA_OBJS)

# The test suite executable recipe
unit_tests: $(PFUNIT_OBJS) $(PFUNIT_DRIVER) $(STELLA_OBJS)
	@echo "Compiling: pFUnit tests"
	$(FC) -o $@ $^ -DUSE_MPI=YES $(PFUNIT_INC) $(PFUNIT_LIB) -I. $(SIMPLEDATAIO_LIB_ABS) $(LIBS) $(F90FLAGS) -I$(PFUNIT_TEST_DIR) -D_TEST_SUITES='"$(PFUNIT_TEST_LIST)"'

check-unit: check-pfunit-install unit_tests
	@echo "Running: pFUnit tests"
	@mpirun -np $(PFUNIT_NPROC) $(PFUNIT_TEST_PATH)/unit_tests

###############################################################
# Coverage

.PHONY: coverage pfunit-coverage run-gcov

coverage:
	gcov *.f90

ifndef USE_GCOV
coverage-pfunit: coverage_gcov_fail
coverage_gcov_fail:
	$(error Cannot make coverage-pfunit as USE_GCOV undefined. Rebuild everything with USE_GCOV=on)
else
coverage-pfunit: stella unit_tests coverage
	mkdir -p coverage-pfunit
	mv *.gcov coverage-pfunit/. || echo "Couldn't move generated gcov files to coverage-pfunit. This is likely because gcov files not found in ${GK_HEAD_DIR}."
endif

###############################################################
# Clean up

.PHONY: clean-unit

clean-unit:
	@echo "Cleaning pFUnit tests"
	$(RM) unit_tests $(PFUNIT_TEST_DIR)/*.o $(PFUNIT_TEST_DIR)/*.mod \
      $(PFUNIT_TEST_DIR)/*~ $(PFUNIT_TEST_DIR)/*.gcov $(PFUNIT_TEST_DIR)/*.gcda \
      $(PFUNIT_TEST_DIR)/*.gcno $(PFUNIT_TEST_DIR)/*.log $(PFUNIT_TEST_DIR)/*.error \
      $(PFUNIT_OBJS) $(PFUNIT_INTERMEDIATES) ./test_*mod ./wraptest_*mod

clean: clean-unit
