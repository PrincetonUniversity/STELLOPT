!-----------------------------------------------------------------------
!     Module:        nescoil_input_mod
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          03/20/2026
!     Description:   This module contains the NESCOIL input namelist and
!                    subroutine which initializes and reads the
!                    NESCOIL input namelist.
!-----------------------------------------------------------------------
      MODULE nescoil_input_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE safe_open_mod, ONLY: safe_open
      USE nescoil_globals, ONLY: w_psurf, w_csurf, w_bnuv, w_jsurf, &
         w_xerr, w_svd, mstrt, mstep, mkeep, mdspw, curwt, trgwt, &
         nu, nv, nu1, nv1, npol, ntor, mf, nf, md, nd, &
         nmax, mnd, nuv, nuv1, nuvh, nuvh1, np, &
         cr, cz, ms, ns, &
         ms1, ns1, cr1, cz1, cl1, iota_edge, &
         phip_edge, curpol, &
         cut, cup, ibex, &
         cr2, cz2, cr3, cz3, cf, sf, inesc
      !USE outctrl, ONLY: w_psurf, w_csurf, w_bnuv, w_jsurf, w_xerr, w_svd
      !use SvdCtrl, ONLY: mstrt, mstep, mkeep, mdspw, curwt, trgwt
      !USE vmeshes, ONLY: nu, nv, nu1, nv1, npol, ntor, mf, nf, md, nd, &
      !                   nmax, mnd, nuv, nuv1, nuvh, nuvh1
      !USE Vprecal1, ONLY: np
      !USE Vvacuum1, ONLY: cr, cz, ms, ns
      !USE Vvacuum2, ONLY: ms1, ns1, cr1, cz1, cl1, iota_edge, &
      !                    phip_edge, curpol
      !USE Vvacuum3, ONLY: cut, cup, ibex
      !USE vvacuum4, ONLY: cr2, cz2 ! not used
      !USE vvacuum5, ONLY: cr3, cz3 ! not used
      !USE vvacuum6, ONLY: cf, sf   ! not used

!-----------------------------------------------------------------------
!     Module Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, PARAMETER :: MMAX_IN = 24
      INTEGER, PARAMETER :: NMAX_IN = 24
      REAL(rprec), DIMENSION(-NMAX_IN:NMAX_IN,0:MMAX_IN) :: &
                                                RBC_SURF, ZBS_SURF, &
                                                RBC_PLASMA, ZBS_PLASMA,&
                                                LBS_PLASMA

!-----------------------------------------------------------------------
!     Input Namelists
!         &nescoil_input
!            nu          Number of potential surface poloidal gridpoints
!            nv          Number of potential surface toroidal gridpoints
!            nu1         Number of equilibrium surface poloidal gridpoints
!            nv1         Number of equilibrium surface toroidal gridpoints
!            npol        Not explicitly used (maybe in accuracy)
!            ntor        Not explicitly used (maybe in accuracy)
!            mf          Poloidal modes in current potential
!            nf          Toroidal modes in current potential
!            md          Poloidal modes in eq/pot surface representation
!            nd          Toroidal modes in eq/pot surface representation
!            np          Number of field periods
!            iota_edge   Edge rotational transform (for min_xerr_phimn)
!            phip_edge   Toroidal flux deriviative (for min_xerr_phimn)
!            curpol      Total poloidal current per field period [A] (not used)
!            cut         Toroidal current flag (toroidal coils)
!            cup         Poloidal current flag (modular coils)
!            ibex        External Current flag (1/R field)
!            mstrt       Method + svdscan start if >1
!                        >=0:Berr, <0:Xerr, unless MSTEP <=0: Least square
!            mstep       Method + svdscan stepsize:
!                        <=0: LeastSquare, =0: use old f04abe (now SOLVER), no svd
!            mkeep       svd/scan control: =0: svdscan, else keep |nkeep| wgts
!                        write all weights to output
!            mdspw       2+exponent of dsur multiplying bfn,ben:
!                        <0: post-calculate Xerr svdscan and write to output
!            curwt       Weight for surface current minimization
!                        Works ONLY in LSQ branch
!            trgwt       Not implemented yet (PMV)
!                 Note for w_ vars (-) implies one, (+) implies all output
!            w_psurf     Write plasma surface info 
!                        1: R/Z, 2: X/Y/Z, 3: NX/NY/NZ, 4: dXdu/dYdu/dXdv/dYdv
!            w_csurf     Write current surface info 
!                        1: R/Z, 2: X/Y/Z, 3: NX/NY/NZ, 4: JX/JY/JZ
!            w_bnuv      Write Bnorm field info
!                        1: See accuracy.f, 2: BN_EXT
!            w_jsurf     Write potential info
!                        1: Potential, 2: Current
!            w_xerr      Write X error (displacement) info (not implemented)
!            w_svd       Write SVD info
!                        2: Weights
!            RBC_PLASMA  Plasma Rmn_cos surface harmonics
!            ZBS_PLASMA  Plasma Zmn_sin surface harmonics
!            LBS_PLASMA  Plasma Lmn_sin (lambda) surface harmonics
!            RBC_SURF    Potential Rmn_cos surface harmonics
!            ZBS_SURF    Potential Zmn_sin surface harmonics
!            
!-----------------------------------------------------------------------
      NAMELIST /nescoil_input/ nu, nv, nu1, nv1, npol, ntor, &
                               mf, nf, md, nd, &
                               np, iota_edge, phip_edge, curpol, &
                               cut, cup, ibex, &
                               mstrt,mstep,mkeep,mdspw,curwt,trgwt, &
                               w_psurf, w_csurf, w_bnuv, w_jsurf, &
                               w_xerr, w_svd, &
                               RBC_PLASMA, ZBS_PLASMA, LBS_PLASMA, &
                               RBC_SURF, ZBS_SURF
      
!-----------------------------------------------------------------------
!     Subroutines
!         init_nescoil_input:   Initializes the namelist
!         read_beams3d_input:   Reads beams3d_input namelist
!-----------------------------------------------------------------------
      CONTAINS

      SUBROUTINE init_nescoil_input
      IMPLICIT NONE
      nu     = -1
      nv     = -1
      nu1     = -1
      nv1     = -1
      npol     = 1
      ntor    = 1
      mf     = -1
      nf     = -1
      md     = -1
      nd     = -1
      np     = 1
      iota_edge = 0.0
      phip_edge = 0.0
      curpol    = 0.0
      cut       = 0
      cup       = 1
      ibex      = 0
      mstrt     = 0
      mstep     = 0
      mkeep     = 0
      mdspw     = 4
      curwt     = 0.0
      trgwt     = 0.0
      w_psurf   = 0
      w_csurf   = 0
      w_bnuv    = 0
      w_jsurf   = 0
      w_xerr    = 0
      w_svd     = 0

      RBC_PLASMA = 0.0
      ZBS_PLASMA = 0.0
      LBS_PLASMA = 0.0
      RBC_SURF   = 0.0
      ZBS_SURF   = 0.0

      RETURN
      END SUBROUTINE init_nescoil_input
      
      SUBROUTINE read_nescoil_input(filename, istat,loutput)
      IMPLICIT NONE
      CHARACTER(*), INTENT(in) :: filename
      INTEGER, INTENT(out) :: istat
      LOGICAL, INTENT(in) :: loutput ! turns output on
      LOGICAL :: lexist
      INTEGER :: iunit, m, n, ntotal
      CHARACTER(LEN=1000) :: line

      ! Read namelist
      IF (filename /= 'IMAS') THEN
         istat=0
         iunit=12
         INQUIRE(FILE=TRIM(filename),EXIST=lexist)
         IF (.not.lexist) THEN
            istat = -1
            RETURN
         END IF
         CALL safe_open(iunit,istat,TRIM(filename),'old','formatted')
         IF (istat /= 0) THEN
            WRITE(6,'(A)') 'ERROR opening file: ',TRIM(filename)
            CALL FLUSH(6)
            RETURN
         END IF
         READ(iunit,NML=nescoil_input,IOSTAT=istat)
         IF (istat /= 0) THEN
            WRITE(6,'(A)') 'ERROR reading namelist NESCOIL_INPUT from file: ',TRIM(filename)
            backspace(iunit)
            read(iunit,fmt='(A)') line
            write(6,'(A)') 'Invalid line in namelist: '//TRIM(line)
            CALL FLUSH(6)
            CLOSE(iunit)
            RETURN
         END IF
         CLOSE(iunit)
      END IF

      IF (loutput) THEN
         CALL safe_open(inesc, istat, 'nescout.'//TRIM(filename(7:)), 'unknown', 'formatted')
         IF (istat .ne. 0) stop 'Error opening nescout file'
      END IF

      ! Now handle NESCOIL variables
      nmax  = 3*nv*(nu+2)+1
      nuv   = nu*nv
      nuv1  = nu1*nv1
      nuvh  = nuv/2 + nu
      nuvh1 = nuv1/2 + nu1

      IF (loutput) THEN
         write(inesc, '(A)') '----- Grid Spatial Dimensions -----'
         write(inesc, '(A)') 'nu, nv, nu1, nv1, npol, ntor'
         write(inesc,"(6i6,l)")  nu, nv, nu1, nv1, npol, ntor, .FALSE.
      END IF

      ! Fix mf and nf
      !md = 0; nd = 0
      DO n = -NMAX_IN, NMAX_IN
         DO m = 0, MMAX_IN
            IF ((RBC_SURF(n,m) /= 0) .or. &
                (ZBS_SURF(n,m) /= 0) .or. &
                (RBC_PLASMA(n,m) /= 0) .or. &
                (ZBS_PLASMA(n,m) /= 0) .or. &
                (LBS_PLASMA(n,m) /= 0)) THEN
               md = MAX(m,md)
               nd = MAX(ABS(n),nd)
            END IF
         END DO
      END DO
      md = MAX(md,2) ! Not sure why this is needed but it is
      nd = MAX(nd,2) ! Not sure why this is needed but it is
      mnd   = (md + 1)*(2*nd + 1)

      IF (loutput) THEN
         write (inesc, '(A)') '----- Fourier Dimensions -----'
         write (inesc, '(A)') 'mf, nf, md, nd'
         write (inesc,"(4i6)")  mf, nf, md, nd
         write (inesc, '(A)') '----- Plasma information from VMEC -----'
         write (inesc, '(A)') 'np, iota_edge, phip_edge, curpol'
         write (inesc,"(i6,3g25.16)")  np, iota_edge, phip_edge, curpol
         write (inesc, '(A)') '----- Current Controls -----'
         write (inesc, '(A)') 'cut, cup, ibex'
         write (inesc,"(2g25.16,i6)")  cut, cup, ibex
         write (inesc, '(A)') '----- SVD controls -----'
         write (inesc, '(A)') 'mstrt, mstep, mkeep, mdspw, curwt, trgwt'
         write (inesc,"(4i6,2g25.16)")  mstrt,mstep,mkeep,mdspw,curwt,trgwt
         write (inesc, '(A)') '----- Output controls -----'
         write (inesc, '(A)') 'w_psurf, w_csurf, w_bnuv, w_jsurf, w_xerr, w_svd'
         write (inesc,"(6i6)") w_psurf,w_csurf,w_bnuv,w_jsurf,w_xerr,w_svd
      END IF

      ! Allocate fourier arrays
      IF (ALLOCATED(cr)) DEALLOCATE(cr)
      IF (ALLOCATED(cz)) DEALLOCATE(cz)
      IF (ALLOCATED(cr1)) DEALLOCATE(cr1)
      IF (ALLOCATED(cz1)) DEALLOCATE(cz1)
      IF (ALLOCATED(cl1)) DEALLOCATE(cl1)
      IF (ALLOCATED(cr2)) DEALLOCATE(cr2)
      IF (ALLOCATED(cz2)) DEALLOCATE(cz2)
      IF (ALLOCATED(cr3)) DEALLOCATE(cr3)
      IF (ALLOCATED(cz3)) DEALLOCATE(cz3)
      IF (ALLOCATED(cf)) DEALLOCATE(cf)
      IF (ALLOCATED(sf)) DEALLOCATE(sf)
      ALLOCATE(cr(0:md,-nd:nd),cz(0:md,-nd:nd))
      ALLOCATE(cr1(0:md,-nd:nd),cz1(0:md,-nd:nd),cl1(0:md,-nd:nd))
      ALLOCATE(cr2(0:md,-nd:nd),cz2(0:md,-nd:nd))
      ALLOCATE(cr3(0:md,-nd:nd),cz3(0:md,-nd:nd))
      ALLOCATE(cf(0:md,-nd:nd),sf(0:md,-nd:nd))

      ! Setup boundary array
      cr = 0.0; cz = 0.0
      ms = 0; ns = 0
      DO n = -nd,nd
         DO m = 0,md
            cr(m,n) = RBC_SURF(n,m)
            cz(m,n) = ZBS_SURF(n,m)
            ms = max(ms,abs(m))
            ns = max(ns,abs(n))
         END DO
      END DO
      cr(0:ms,0) = .5_rprec*cr(0:ms,0)
      cz(0:ms,0) = .5_rprec*cz(0:ms,0)

      ! Setup plasma array
      cr1 = 0.0; cz1 = 0.0; cl1 = 0.0
      ms1 = 0; ns1 = 0
      DO n = -nd,nd
         DO m = 0,md
            cr1(m,n) = RBC_PLASMA(n,m)
            cz1(m,n) = ZBS_PLASMA(n,m)
            cl1(m,n) = LBS_PLASMA(n,m)
            ms1 = max(ms1,abs(m))
            ns1 = max(ns1,abs(n))
         END DO
      END DO
      cr1(0:ms1,0) = .5_rprec*cr1(0:ms1,0)
      cz1(0:ms1,0) = .5_rprec*cz1(0:ms1,0)

      IF (loutput) THEN
         ntotal = COUNT((RBC_PLASMA /= 0) .or. (ZBS_PLASMA /=0) .or. (LBS_PLASMA /=0))
         write (inesc, '(A)') '----- Plasma Surface -----'
         write (inesc, '(A)') 'Number of fourier modes in table'
         write (inesc, *)  ntotal
         write (inesc, '(A)') '----- Plasma boundary fourier coefficients  -----'
         write (inesc, '(A)') '   m    n        R(m,n)     Z(m,n)    Lamda(m,n)'
         DO n = -NMAX_IN, NMAX_IN
            DO m = 0, MMAX_IN
               IF ((RBC_PLASMA(n,m) /= 0) .or. &
                   (ZBS_PLASMA(n,m) /= 0) .or. &
                   (LBS_PLASMA(n,m) /= 0)) THEN
                  write (inesc,"(2i4,6g20.10)") &
                     m, n, RBC_PLASMA(n,m), ZBS_PLASMA(n,m), LBS_PLASMA(n,m),0.0,0.0,0.0
               END IF
            END DO
         END DO
         ntotal = COUNT((RBC_SURF /= 0) .or. (ZBS_SURF/=0))
         write (inesc, '(A)') '----- Coil Surface -----'
         write (inesc, '(A)') 'Number of fourier modes in table'
         write (inesc, *)  ntotal
         write (inesc, '(A)') '----- Coil surface fourier coefficients -----'
         write (inesc, '(A)') '    m    n         R(m,n)         Z(m,n)'
         DO n = -NMAX_IN, NMAX_IN
            DO m = 0, MMAX_IN
               IF ((RBC_SURF(n,m) /= 0) .or. &
                   (ZBS_SURF(n,m) /= 0)) THEN
                  write (inesc,"(2i4,4g20.10)") m, n, RBC_SURF(n,m), ZBS_SURF(n,m), 0.0, 0.0
               END IF
            END DO
         END DO
         write (inesc, '(A)') '----- end inputs, begin outputs. Nescoil Version 1.0 -----'
         if( (ABS(cup) == 0) .and. (ABS(cut) == 0) ) then
            write (inesc, '(A)') '----- Solving for Saddle coils -----'
         else
            write (inesc, '(A)') '----- Solving for Modular coils -----'
         endif

         if( ibex .eq. 0 ) then
            write (inesc, '(A)') '----- No background coils used -----'
         else
            write (inesc, '(A)') '----- Background coils used -----'
         endif
      END IF

      RETURN

      END SUBROUTINE read_nescoil_input

      SUBROUTINE write_nescoil_namelist(iunit_out, istat)
      INTEGER, INTENT(in) :: iunit_out
      INTEGER, INTENT(out) :: istat
      INTEGER :: n, m
      CHARACTER(LEN=*), PARAMETER :: outboo  = "(2X,A,1X,'=',1X,L1)"
      CHARACTER(LEN=*), PARAMETER :: outint  = "(2X,A,1X,'=',1X,I0)"
      CHARACTER(LEN=*), PARAMETER :: outflt  = "(2X,A,1X,'=',1X,ES22.12E3)"
      CHARACTER(LEN=*), PARAMETER :: outexp  = "(2X,A,1X,'=',1X,ES22.12E3)"
      CHARACTER(LEN=*), PARAMETER :: outcmp  = "(2x,A,1X,'=','(',i3,',',i3,')')"
      CHARACTER(LEN=*), PARAMETER :: outstr  = "(2X,A,1X,'=',1X,'''',A,'''')"
      CHARACTER(LEN=*), PARAMETER :: onevar  = "(2X,A,1X,'=',1X,L1,2(2X,A,1X,'=',1X,ES22.12E3))"
      CHARACTER(LEN=*), PARAMETER :: vecvar  = "(2X,A,'(',I3.3,')',1X,'=',1X,ES22.12E3)"
      CHARACTER(LEN=*), PARAMETER :: vecvar2  = "(2X,A,'(',I3.3,',',I3.3,')',1X,'=',1X,ES22.12E3)"
      istat = 0
      WRITE(iunit_out,'(A)') '&NESCOIL_INPUT'
      WRITE(iunit_out,'(A)') '!---------- Plasma Surface Parameters ------------'
      WRITE(iunit_out,outint) 'NU1',nu1
      WRITE(iunit_out,outint) 'NV1',nv1
      WRITE(iunit_out,outint) 'MD',md
      WRITE(iunit_out,outint) 'ND',nd
      WRITE(iunit_out,outint) 'NP',np
      WRITE(iunit_out,outflt) 'IOTA_EDGE',iota_edge
      WRITE(iunit_out,outflt) 'PHIP_EDGE',phip_edge
      WRITE(iunit_out,outflt) 'CURPOL',curpol
      WRITE(iunit_out,'(A)') '!---------- Potential Surface Parameters ------------'
      WRITE(iunit_out,outint) 'NU',nu
      WRITE(iunit_out,outint) 'NV',nv
      WRITE(iunit_out,outint) 'MF',mf
      WRITE(iunit_out,outint) 'NF',nf
      WRITE(iunit_out,outflt) 'CUT',cut
      WRITE(iunit_out,outflt) 'CUP',cup
      WRITE(iunit_out,outint) 'IBEX',ibex
      WRITE(iunit_out,'(A)') '!---------- Solver Parameters ------------'
      WRITE(iunit_out,outint) 'MSTRT',mstrt
      WRITE(iunit_out,outint) 'MSTEP',mstep
      WRITE(iunit_out,outint) 'MKEEP',mkeep
      WRITE(iunit_out,outint) 'MDSPW',mdspw
      WRITE(iunit_out,outflt) 'CURWT',curwt
      WRITE(iunit_out,outflt) 'TRGWT',trgwt
      WRITE(iunit_out,'(A)') '!---------- Output Parameters ------------'
      WRITE(iunit_out,outint) 'W_PSURF',w_psurf
      WRITE(iunit_out,outint) 'W_CSURF',w_csurf
      WRITE(iunit_out,outint) 'W_BNUV',w_bnuv
      WRITE(iunit_out,outint) 'W_JSURF',w_jsurf
      WRITE(iunit_out,outint) 'W_XERR',w_xerr
      WRITE(iunit_out,outint) 'W_SVD',w_svd
      WRITE(iunit_out,'(A)') '!---------- Equilibrium Surface Harmonics ------------'
      DO m = 0, md
         DO n = -nd, nd
            IF ((cr1(m,n).ne.0) .or. (cz1(m,n).ne.0) .or. (cl1(m,n).ne.0)) THEN
               WRITE(iunit_out,'(3(A,I4.3,A,I3.3,A,ES22.12E3))') &
                  '  RBC_PLASMA(',n,',',m,') = ',cr1(m,n), &
                  '    ZBS_PLASMA(',n,',',m,') = ',cz1(m,n), &
                  '    LBS_PLASMA(',n,',',m,') = ',cl1(m,n)
            END IF
         END DO
      END DO
      WRITE(iunit_out,'(A)') '!---------- Potential Surface Harmonics ------------'
      DO m = 0, md
         DO n = -nd, nd
            IF ((cr(m,n).ne.0) .or. (cz(m,n).ne.0)) THEN
               WRITE(iunit_out,'(2(A,I4.3,A,I3.3,A,ES22.12E3))') &
                  '  RBC_SURF(',n,',',m,') = ',cr(m,n), &
                  '    ZBS_SURF(',n,',',m,') = ',cz(m,n)
            END IF
         END DO
      END DO
      WRITE(iunit_out,'(A)') '/'

      END SUBROUTINE write_nescoil_namelist

      SUBROUTINE write_nescoil_namelist_byfile(filename)
      CHARACTER(LEN=*), INTENT(in) :: filename
      INTEGER :: iunit, istat
      LOGICAL :: lexists
      
      iunit = 100
      istat = 0
      INQUIRE(FILE=TRIM(filename),exist=lexists)
      IF (lexists) THEN
         OPEN(unit=iunit, file=TRIM(filename), iostat=istat, status="old", position="append")
      ELSE
         OPEN(unit=iunit, file=TRIM(filename), iostat=istat, status="new")
      END IF
      IF (istat .ne. 0) RETURN
      CALL write_nescoil_namelist(iunit,istat)
      CLOSE(iunit)

      RETURN
      END SUBROUTINE write_nescoil_namelist_byfile

      END MODULE nescoil_input_mod
