!Main program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM KNOSOS

  USE GLOBAL
  USE KNOSOS_STELLOPT_MOD
#ifdef MPIandPETSc
  USE petscsys
#endif
  IMPLICIT NONE
  !Parameters
  INTEGER, PARAMETER :: ntime=100000
  !Allocatable
  INTEGER, ALLOCATABLE :: rank(:,:)
  REAL*8, ALLOCATABLE :: nb(:,:,:),dnbdpsi(:,:,:),Tb(:,:,:),dTbdpsi(:,:,:),Epsi(:,:)
  REAL*8, ALLOCATABLE :: Gb(:,:,:),Qb(:,:,:),Sb(:,:,:),Pb(:,:,:)
  !Other
  CHARACTER*100 serr
  INTEGER itime,is,ns,nbb,regb(nbx),jerr
  REAL*8 dt,S(nsx),Zb(nbx),Ab(nbx),fracb(nbx)
#ifdef MPIandPETSc
  INTEGER ierr
#endif
#ifdef IPPoNIFS
! INCLUDE "mpif.h"
#include <petsc/finclude/petscsys.h>
#endif

  !Time
  CHARACTER*30, PARAMETER :: routine="KNOSOS"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)
  !Initialize MPI
  CALL INITIALIZE_MPI()
  !Read input files (simulation parameters, models, flux-surfaces, species...)
  KN_STELLOPT=.FALSE.
  KN_EXT=""
  KN_DKESFILE="boozmn.nc"
  CALL READ_INPUT(dt,ns,s,nbb,Zb,Ab,regb,fracb)
  !Allocate some transport-related quantities
  ALLOCATE(nb(nbb,ns,nerr),dnbdpsi(nbb,ns,nerr),Tb(nbb,ns,nerr),dTbdpsi(nbb,ns,nerr),&
    & Epsi(ns,nerr),Gb(nbb,ns,nerr),Qb(nbb,ns,nerr),Sb(nbb,ns,nerr),Pb(nbb,ns,nerr),rank(ns,nerr))
  !Initialize PETSC and the random number generator according to the input, distribute MPI jobs
  CALL DISTRIBUTE_MPI(ns,rank)
#ifdef MPIandPETSc
  CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,myrank,myrank,PETSC_COMM_WORLD,ierr)
  CALL PETSCINITIALIZE(PETSC_NULL_CHARACTER,ierr)
#endif
  IF(nerr.GT.0.OR.FAST_IONS) CALL INIT_RANDOMSEED(0)
  WRITE(iout,*) 'Initialized'
  CALL INIT_FILES()
  !Perform a loop in samples
  DO jerr=1,nerr
     !Follow time evolution  
     DO itime=1,ntime
        !Perform a loop in flux-surfaces
        DO is=1,ns
           !Êach process has been asigned some particular case
           IF(myrank.NE.rank(is,jerr)) CYCLE
           !Read the magnetic field
           IF(itime.EQ.1) THEN
              CALL READ_BFIELD(s(is))
              !Create (or read) a DKES-like database of monoenergetic transport coefficients
              CALL CALC_DATABASE(s,is,ns)
              IF(KNOSOS_STELLOPT) WRITE(600+myrank,'(12(1pe13.5))') s(is),&
                & KN_1NU,KN_SNU,KN_SBP,KN_GMC,KN_GMA,KN_QER,KN_VB0,KN_VBB,KN_WBW,KN_DBO,KN_VBM,KN_FTR
              IF(ONLY_DB) CYCLE
              !Read the plasma profiles 
              CALL READ_PLASMAS(nbb,fracb(1:nbb),s(is),Zb(1:nbb),Ab(1:nbb),nb(:,is,jerr),dnbdpsi(:,is,jerr),&
                   & Tb(:,is,jerr),dTbdpsi(:,is,jerr),Epsi(is,jerr)) 
              !Read the particle and energy sources
              CALL READ_SOURCES(nbb,s(is),Sb(:,is,jerr),Pb(:,is,jerr))
           END IF
           !Calculate the flux for each surface
           !(incl. drift-kinetic eq., ambipolarity and quasineutrality)
           CALL SOLVE_DKE_QN_AMB(itime,nbb,Zb(1:nbb),Ab(1:nbb),regb(1:nbb),s(is),&
                & nb(:,is,jerr),dnbdpsi(:,is,jerr),Tb(:,is,jerr),dTbdpsi(:,is,jerr),&
                & Epsi(is,jerr),Gb(:,is,jerr),Qb(:,is,jerr))
        END DO
        !Evolve the plasma profiles (without error propagation)
        CALL TRANSPORT(dt,nbb,ns,s(1:ns),Zb(1:nbb),Ab(1:nbb),REGB(1:nbb),&
             & nb(:,:,jerr),dnbdpsi(:,:,jerr),Gb(:,:,jerr),Sb,&
             & Tb(:,:,jerr),dTbdpsi(:,:,jerr),Qb(:,:,jerr),Pb,Epsi(:,jerr))
        IF(dt.LT.0.AND..NOT.SS_IMP) EXIT
     END DO
  END DO
  IF(nerr.GT.1) CALL AVERAGE_SAMPLES(nbb,ns,s(1:ns),Epsi,Gb,Qb) 
  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)
  !End MPI and PETSC
  serr="Simulation complete"
  CALL END_ALL(serr,.TRUE.)

END PROGRAM KNOSOS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


