       SUBROUTINE summodosd(nsurf, npt, x, c2, k_s, bsuppar,
     1            b2, iotac, iotapc, prespc)
       USE stel_kinds
       USE normalize_data, ONLY: lasym_v                                      ! 110909 RS: Logical for asymmetric input (If true)
       USE ballooning_data
       USE general_dimensions
       USE fmesh_quantities
       USE summod
       IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       INTEGER,INTENT(in) :: nsurf, npt
       REAL(rprec),INTENT(in), DIMENSION (npt) :: x
       REAL(rprec),INTENT(in):: iotac, prespc, iotapc
       REAL(rprec),INTENT(out), DIMENSION(npt) :: c2,
     1   k_s, b2, bsuppar
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
       REAL(rprec), PARAMETER :: one = 1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER :: j, lj
       REAL(rprec) :: zetacn, thetacn, twopi, alpha,
     1  zetacn_i, rdum
!------------------------------------------------------------------------

       CALL alloc_summod(npt)
       twopi   = 8*ATAN(one)

       IF (l_geom_input) THEN                                                ! Location of initial field line given in (THETA, ZETA) format. Convert to radians.
 
         zetacn  = (twopi*init_zeta)/360._dp                                 ! Convert init_zeta [= zeta_k] to radians; notice it CANNOT exceed 2*pi
         thetacn = (twopi*init_theta)/360._dp                                ! Convert init_theta to radians....
         CALL obtain_field_line(zetacn, thetacn, iotac, nsurf,               ! ... and find label (ALPHA) of line that contains (init_theta, init_zeta)
     1     alpha)

       ELSE                                                                  ! Location of initial field line given in (ALPHA, ANGLE) format.

         IF (l_tokamak_input) THEN                                           ! TOKAMAK-like input => ALPHA=ALPHA_TOK; ANGLE=THETA. COBRA converts to stellarator params
                                                                             !   Beware that this might be problematic, since phi(theta) may not be uni-valued in VMEC coordinates
           alpha = iotac*(twopi*init_alpha_tok)/360._dp                      !   Standard ALPHA_TOK = Q*THETA +Q*LAMBDA - ZETA = ALPHA_ST/IOTA; Convert to stellarator ALPHA in radians
           thetacn = (twopi*init_thetak_tok)/360._dp                         ! Convert int_thetak_tok [=Tokamak theta_k] to radians;  notice it CANNOT exceed 2*pi
           zetacn_i = (thetacn - alpha)/iotac                                ! Use as guess to find real zeta_k..
           CALL obtain_zeta(alpha, iotac, zetacn_i, nsurf, 1,                ! ...using OBTAIN_ZETA
     1       zetacn, thetacn, rdum)

         ELSE                                                                ! STELLARATOR-like input => ALPHA=ALPHA_ST; ANGLE=ZETA

           alpha = (twopi*init_alpha_st)/360._dp                             ! Convert alpha to radians; notice it CANNOT exceed 2*pi
           zetacn  = (twopi*init_zetak_st)/360._dp                           ! Convert init_zetak_st [= zeta_k] to radians; notice it CANNOT exceed 2*pi
           thetacn = alpha + iotac*zetacn                                    ! Use this as first guess for OBTAIN_THETA subroutine

         ENDIF

       ENDIF

       zetang  = zetacn+x                                                    ! Recall X = ZETA-ZETA_k; P, Q and R MUST HOWEVER be evaluated at ZETA = X+ZETA_k
       CALL obtain_theta(alpha, iotac, thetacn, nsurf, npt,                  ! Once in stellarator (ALPHA, ZETA) format, compute THETAs for all ZETA's along line used in calculation.
     1   zetang, thetang, lambdath)

       IF (lfail_balloon) THEN
          CALL free_summod
          RETURN
       END IF
!=============
!   BEGIN FOURIER INVERSION
!=============

!...   initialize for Fourier inversion

       bfield = 0; bfieldze = 0; bfields = 0; bfieldth = 0
       rboo = 0; rs = 0; rze = 0; rth = 0; zs = 0; zze = 0
       zth = 0; lambdaze = 0; lambdas = 0; bsupze = 0; bsupth = 0

       fourier1: DO j = 1, mnmax_v                                            ! Fourier invert back to real space

         arg = xm_v(j)*thetang-zetang*xn_v(j)
         ccosi = COS(arg)
         ssine = SIN(arg)
         lj = mnmax_v*(nsurf-1)+j
         rboo = rboo+rmncf(lj)*ccosi                                         ! cylindrical R
         IF (lasym_v) rboo = rboo+rmnsf(lj)*ssine                            ! 110909 RS: Asymmetric input
         rth = rth-rmncf(lj)*xm_v(j)*ssine                                   ! ..... theta derivative
         IF (lasym_v) rth = rth+rmnsf(lj)*xm_v(j)*ccosi                      ! 110909 RS: Asymmetric input
         rze = rze+rmncf(lj)*xn_v(j)*ssine                                   ! ..... zeta derivative
         IF (lasym_v) rze = rze-rmnsf(lj)*xn_v(j)*ccosi                      ! 110909 RS: Asymmetric input
         rs = rs+rmncpf(lj)*ccosi                                            ! ..... radial derivative
         IF (lasym_v) rs = rs+rmnspf(lj)*ssine                               ! 110909 RS: Asymmetric input
         zth = zth+zmnsf(lj)*xm_v(j)*ccosi                                   ! cylindrical Z: theta derivative
         IF (lasym_v) zth = zth-zmncf(lj)*xm_v(j)*ssine                      ! 110909 RS: Asymmetric input
         zze = zze-zmnsf(lj)*xn_v(j)*ccosi                                   ! ..... zeta derivative
         IF (lasym_v) zze = zze+zmncf(lj)*xn_v(j)*ssine                      ! 110909 RS: Asymmetric input
         zs = zs+zmnspf(lj)*ssine                                            ! ..... radial derivative
         IF (lasym_v) zs = zs+zmncpf(lj)*ccosi                               ! 110909 RS: Asymmetric input
         lambdas = lambdas+lmnspf(lj)*ssine                                  ! lambda radial derivative
         IF (lasym_v) lambdas = lambdas+lmncpf(lj)*ccosi                     ! 110909 RS: Asymmetric input
         lambdaze = lambdaze-xn_v(j)*lmnsf(lj)*ccosi                         ! ..... zeta derivative
         IF (lasym_v) lambdaze = lambdaze+xn_v(j)*lmncf(lj)*ssine            ! 110909 RS: Asymmetric input

       ENDDO fourier1

       fourier2: DO j = 1, mnmax_vnyq                                            ! Fourier invert back to real space

         arg = xm_vnyq(j)*thetang-zetang*xn_vnyq(j)
         ccosi = COS(arg)
         ssine = SIN(arg)
         lj = mnmax_vnyq*(nsurf-1)+j
         bfield = bfield+bmncf(lj)*ccosi                                     !magnetic field magnitude
         IF (lasym_v) bfield = bfield + bmnsf(lj)*ssine                      ! 110909 RS : Asymmetric input
         bfields = bfields+bmncpf(lj)*ccosi                                  ! ..... radial derivative
         IF (lasym_v) bfields = bfields + bmnspf(lj)*ssine                   ! 110909 RS : Asymmetric input
         bfieldze = bfieldze+xn_vnyq(j)*bmncf(lj)*ssine                         ! ..... zeta derivative
         IF (lasym_v) bfieldze = bfieldze-xn_vnyq(j)*bmnsf(lj)*ccosi            ! 110909 RS: Asymmetric input
         bfieldth = bfieldth-xm_vnyq(j)*bmncf(lj)*ssine                         ! ..... theta derivative
         IF (lasym_v) bfieldth = bfieldth+xm_vnyq(j)*bmnsf(lj)*ccosi            ! 110909 RS: Asymmetric input
         bsupth= bsupth+ bsupumncf(lj)*ccosi                                 ! contravariant theta-comp. magnetic field
         IF (lasym_v) bsupth= bsupth+ bsupumnsf(lj)*ssine                    ! 110909 RS: Asymmetric input
         bsupze= bsupze+ bsupvmncf(lj)*ccosi                                 ! contravariant zeta-comp. magnetic field
         IF (lasym_v) bsupze= bsupze+ bsupvmnsf(lj)*ssine                    ! 110909 RS: Asymmetric input

       ENDDO fourier2

!============
!   END FOURIER INVERSION
!============

!...   auxiliar quantities

       rboo2 = rboo**2
       bfieldi = one/bfield
       b2 = bfield**2
       bfield2i = bfieldi**2

!...   VMEC lower metric elements

       gtssub = rth*rs+zth*zs
       gstsub = gtssub
       gzssub = rze*rs+zze*zs
       gszsub = gzssub
       gtzsub = rth*rze+zth*zze
       gztsub = gtzsub
       gttsub = (rth**2)+(zth**2)
       gsssub = (rs**2)+(zs**2)
       gzzsub = (rze**2)+(zze**2)+rboo2

!...    VMEC jacobian

       jacob2 = gsssub*gttsub*gzzsub +
     1    2*gstsub*gtzsub*gszsub -
     2    gttsub*gszsub**2 -gzzsub*gstsub**2-
     3    gsssub*gztsub**2
       rjac2i = one/jacob2

!...   VMEC upper metric elements

       gsssup = (gttsub*gzzsub-gtzsub*gztsub)*rjac2i
       gttsup = (gsssub*gzzsub-gszsub*gzssub)*rjac2i
       gzzsup = (gsssub*gttsub-gstsub*gtssub)*rjac2i
       gstsup = (gztsub*gzssub-gstsub*gzzsub)*rjac2i
       gtssup = gstsup
       gszsup = (gtssub*gtzsub-gttsub*gzssub)*rjac2i
       gzssup = gszsup
       gtzsup = (gstsub*gzssub-gsssub*gztsub)*rjac2i
       gztsup = gtzsup

!...   covariant B-components

       bsubs=  bsupze*gszsub+bsupth*gstsub
       bsubth=  bsupze*gtzsub+bsupth*gttsub
       bsubze=  bsupze*gzzsub+bsupth*gtzsub

!...   VMEC covariant curvature components

       aux = (bsupze*bfieldze+bsupth*bfieldth)*bfield2i
       cks_vmec = bfieldi*(bfields + bfieldi*prespc - bsubs*aux)
       ckth_vmec = bfieldi*(bfieldth - bsubth*aux)

!...   normal curvature in (s,alpha,phi)-coordinates!

       lam1 = 1 + lambdath
       lam2 = - iotac + lambdaze
       lam3 = - iotapc*x+ lambdas
       k_s = cks_vmec - ckth_vmec*lam3/lam1

!...   normal vector squared

       c2 = gzzsup*lam2**2 + gttsup*lam1**2 +
     1   gsssup*lam3**2 + 2*lam2*lam3*gzssup +
     1   2*lam3*lam1*gtssup + 2*lam1*lam2*gztsup

!...   Contravariant zeta-component of field going to B.grad!!

       bsuppar = bsupze

       CALL free_summod

       END SUBROUTINE summodosd
