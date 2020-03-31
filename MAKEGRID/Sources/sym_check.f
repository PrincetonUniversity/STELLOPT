      MODULE sym_check
      
      USE biotsavart
      USE makegrid_global, ONLY: rprec
      USE write_mgrid, ONLY: rmin, rmax, zmin, zmax
      
      IMPLICIT NONE
      
      PRIVATE

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: sym_nfp, sym_nextcur
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: 
     1 sym_br, sym_bz, sym_bp, sym_b
!-----------------------------------------------
!   P u b l i c  S u b r o u t i n e s  A n d  V a r i a b l e s
!-----------------------------------------------
      PUBLIC :: check_symmetry
      PUBLIC :: init_symmetry
      PUBLIC :: cleanup_symmetry
      
      INTEGER, PUBLIC :: sym_ir, sym_jz, sym_kp
      LOGICAL, PUBLIC :: sym_perform_tests
!-----------------------------------------------

      CONTAINS
      
!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------

      SUBROUTINE init_symmetry
      
      sym_nfp = nfp_bs
      sym_nextcur = SIZE(coil_group)
      
      ALLOCATE( sym_br(sym_ir, sym_jz, sym_kp*sym_nfp) )
      ALLOCATE( sym_bz(sym_ir, sym_jz, sym_kp*sym_nfp) )
      ALLOCATE( sym_bp(sym_ir, sym_jz, sym_kp*sym_nfp) )
      ALLOCATE( sym_b(sym_ir, sym_jz, sym_kp*sym_nfp) )
      
      END SUBROUTINE init_symmetry
      
!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------

      SUBROUTINE cleanup_symmetry
      
      DEALLOCATE(sym_bp)
      DEALLOCATE(sym_bz)
      DEALLOCATE(sym_br)
      DEALLOCATE(sym_b)
      
      END SUBROUTINE cleanup_symmetry

!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------

      SUBROUTINE check_symmetry(ig)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: ig
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, j, ki, fp
      REAL(rprec) :: maxerr, avgerr
!----------------------------------------------- 
 
      IF(sym_perform_tests) THEN
 
            CALL sym_compute_bfield(ig)
            
            WRITE(*,225)
            DO ki=1, sym_kp
                  CALL sym_checkrot(ki,maxerr,avgerr)
                  WRITE(*,300) ki, maxerr, avgerr
            END DO

            WRITE(*,275)
            DO fp=1, sym_nfp
                  CALL sym_checkstel(fp,maxerr,avgerr)
                  WRITE(*,300) fp, maxerr, avgerr
            END DO
      
      END IF
     
 225  FORMAT(/" Error factors, rotational symmetry ",
     1 "(corresponding to KI'th plane):",/,
     2 " KI", 6X, "Max Error", 5X, "Avg Error")
 
 275  FORMAT(/" Error factors, stellarator symmetry ",
     1 "(for FP'th field period):",/,
     2 " FP", 6X, "Max Error", 5X, "Avg Error")
     
 300  FORMAT(I2,4X,ES14.5,ES14.5)
            
      END SUBROUTINE check_symmetry
      
!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------
      SUBROUTINE sym_compute_bfield(ig)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: ig
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: numcoils, i, j, k
      REAL(rprec) :: rrr,ppp,zzz,delr,delp,delz
!----------------------------------------------- 
      
      delr = (rmax - rmin)/(sym_ir - 1)
      delz = (zmax - zmin)/(sym_ir - 1)
      delp = (8*DATAN(1.0D0))/(sym_nfp*sym_kp)
      
      sym_bz = 0D0
      sym_bp = 0D0
      sym_br = 0D0
      
      ppp = 0D0
      DO k=1,sym_kp*sym_nfp
            zzz = zmin
            
            DO j=1,sym_jz
                  rrr = rmin
                  
                  DO i=1,sym_ir
                  ! Compute br, bp, bz at each grid point
                  CALL bfield (rrr, ppp, zzz, sym_br(i,j,k),
     1 sym_bp(i,j,k), sym_bz(i,j,k), ig)
                  
                  ! Compute magnitude, b, at each grid point
                  sym_b(i,j,k) = 
     1 DSQRT(sym_br(i,j,k)**2 + sym_bp(i,j,k)**2 + sym_bz(i,j,k)**2)

                  rrr = rrr + delr
                  END DO
                  zzz = zzz + delz
            END DO
            ppp = ppp + delp
      END DO

      END SUBROUTINE sym_compute_bfield
      
!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------
      
      SUBROUTINE sym_checkrot(ki,maxerr,avgerr)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: ki  ! 1 <= ki <= sym_kp
      REAL(rprec), INTENT(out) :: maxerr, avgerr
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i,j,k,l
      REAL(rprec) :: dbr, dbp, dbz, currerr, avgb
!-----------------------------------------------

      maxerr = 0D0
      avgerr = 0D0
      
      ! Compute maximum and average rotational error factors 
      ! for the set of grids corresponding to k=ki
      DO k=ki+sym_kp, sym_kp*sym_nfp, sym_kp
            currerr = 0D0
            DO j=1, sym_jz
            DO i=1, sym_ir
                        
                  dbr = sym_br(i,j,k) - sym_br(i,j,ki)
                  dbp = sym_bp(i,j,k) - sym_bp(i,j,ki)
                  dbz = sym_bz(i,j,k) - sym_bz(i,j,ki)
                  avgb = (sym_b(i,j,k) + sym_b(i,j,ki))/2.0D0
                        
                  IF(avgb .NE. 0D0) THEN
                  currerr = currerr + DSQRT(dbr**2+dbp**2+dbz**2)/avgb
                  END IF
            END DO
            END DO
                  
            avgerr = avgerr + currerr
            IF(currerr .GT. maxerr) maxerr = currerr
      END DO
      
      IF (sym_nfp .GT. 1) avgerr = avgerr/(sym_nfp-1)
            
      END SUBROUTINE sym_checkrot
      
!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------
      
      SUBROUTINE sym_checkstel(fp,maxerr,avgerr)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: fp  ! 1 <= fp <= sym_nfp
      REAL(rprec), INTENT(out) :: maxerr, avgerr
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i,j,k,h,l, startk, midk, endk
      REAL(rprec) :: dbr, dbp, dbz, avgb, currerr
!-----------------------------------------------

      startk = sym_kp*(fp-1)+1
      midk = startk+sym_kp/2-1
      endk = sym_kp*fp+1
      
      maxerr = 0D0
      avgerr = 0D0
      
      ! Compute maximum stellarator symmetry error factor in the
      ! fp'th field period
      DO k=startk, midk
            l = MOD(endk-(k-startk)-1,sym_kp*sym_nfp)+1
            currerr = 0D0
            DO i=1, sym_ir
            DO j=1, sym_jz
                  h = sym_jz - (j-1)
                  
                  dbr = sym_br(i,j,k)+sym_br(i,h,l)
                  dbp = sym_bp(i,j,k)-sym_bp(i,h,l)
                  dbz = sym_bz(i,j,k)-sym_bz(i,h,l)
                  avgb = (sym_b(i,j,k)+sym_b(i,h,l))/2.0D0
                  
                  IF(avgb .NE. 0D0) THEN
                  currerr = currerr + DSQRT(dbr**2+dbp**2+dbz**2)/avgb
                  END IF
            END DO
            END DO
            
            avgerr = avgerr + currerr
            IF(currerr .GT. maxerr) maxerr = currerr
      END DO
      
      IF(midk-startk+1 .GT. 0) avgerr = avgerr/(midk-startk+1)
      
      END SUBROUTINE sym_checkstel
      
!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------

      END MODULE sym_check
