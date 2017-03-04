      MODULE bivariate
      USE stel_kinds
      USE v3post_rfun, ONLY: ideqfile
      INTEGER :: nrgrid, nzgrid, ndim
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: irz11_bi, irz12_bi
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: irz21_bi, irz22_bi
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: w11_bi, w12_bi
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: w21_bi, w22_bi

      PRIVATE :: irz11_bi, irz12_bi, irz21_bi, irz22_bi,
     &             w11_bi,   w12_bi,   w21_bi,   w22_bi
C-----------------------------------------------
!
!     Written by S. P. Hirshman (01/30/03)
!
!     Local variables:
!
!     irz11_bi(ns*nu,kp) : lower left corner indices on a
!                          rectangular rgrid,zgrid enclosing 
!                          the point (s,u,v)
!     irz22_bi(ns*nu,kp) : upper right corner indices
!     wij_bi(ns*nu,kp)  : weight factors for 4 pt interpolation
!
C-----------------------------------------------
     
      CONTAINS

      SUBROUTINE setbivariate (
     &	rsu, zsu, rgrid, zgrid, kin, kp, ib, ierr)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: kin, kp, ib
      INTEGER, INTENT(out) :: ierr
      REAL(rprec), INTENT(in) :: rsu(:,:), zsu(:,:), rgrid(:), zgrid(:)
      INTEGER :: ku, js, istat, ns, nu, index1d, ir, jz, ir1, jz1
      REAL(rprec) :: rad0, zee0, ri, zj, pr, qz, temp, delr, delz
C-----------------------------------------------
!
!     Compute the indices (ir_bi, jz_bi) and weight factors (wij_bi) 
!     for performing bivariate (4 pt) interpolation FROM a rectangular (R X Z)
!     grid TO a general (non-orthogonal) grid (s, u)
!
!     rgrid(1:nr)     :  r coordinate on a fixed, equally spaced grid 
!     zgrid(1:nz)     :  z coordinate on a fixed, equally spaced grid
!                        (may be different mesh size than rgrid)
!     rsu(1:ns,1:nu)  :  r coordinate at s,u
!     zsu(1:ns,1:nu)  :  z coordinate at s,u
!
!     kin             :  toroidal plane index (1 <= kin <= kp)
!     kp              :  number of toroidal planes
!     ib              :  starting radial s index
!
!     index1d         :  index for making indices/weight arrays 1D (in each toroidal plane)
!                        (improves speed by as much as 2X)
!     
C-----------------------------------------------
      nrgrid = SIZE(rgrid); nzgrid = SIZE(zgrid)
      ns = SIZE(rsu,1);  nu = SIZE(rsu,2)
      ndim = (ns-ib+1)*nu
      istat = 0
      ierr=0

      delr = rgrid(2) - rgrid(1)
      delz = zgrid(2) - zgrid(1)

      IF (.not.ALLOCATED(irz11_bi)) ALLOCATE (irz11_bi(ndim,kp), 
     &    irz12_bi(ndim,kp), irz21_bi(ndim,kp), irz22_bi(ndim,kp),
     &    w11_bi(ndim,kp), w12_bi(ndim,kp), 
     &    w21_bi(ndim,kp), w22_bi(ndim,kp),stat = istat)
      IF (istat .ne. 0) STOP 'Allocation error in setbivariate'

      index1d = 0

      DO ku = 1, nu
      DO js = ib, ns

         index1d = index1d + 1
!
!       CHECK THAT BOUNDARY POINTS ARE INSIDE GRID.  IF NOT, STOP!
!       DO NOT stop in v3post_sub calls; do something else.
         rad0 = rsu(js, ku)
         zee0 = zsu(js, ku)

         IF (rad0.lt.rgrid(1) .or. rad0.gt.rgrid(nrgrid) .or.
     &       zee0.lt.zgrid(1) .or. zee0.gt.zgrid(nzgrid)) THEN
!           PRINT*, 'V3POST-Plasma point outside response function grid!'
!           PRINT*," For : ",TRIM(ideqfile)
           ierr=-1
           RETURN
         ENDIF
!
!       DETERMINE INTEGER INDICES (IR,JZ) FOR LOWER LEFT R, Z CORNER GRID POINT
!
         ir = INT((rad0 - rgrid(1))/delr) + 1
         jz = INT((zee0 - zgrid(1))/delz) + 1
         ir1 = MIN(nrgrid,ir+1)
         jz1 = MIN(nzgrid,jz+1)

!
!        STORE INDICES IN 1D ARRAYS
!
         irz11_bi(index1d,kin) = ir  + nrgrid*(jz-1)
         irz22_bi(index1d,kin) = ir1 + nrgrid*(jz1-1)
         irz12_bi(index1d,kin) = ir  + nrgrid*(jz1-1)
         irz21_bi(index1d,kin) = ir1 + nrgrid*(jz-1)
!
!       COMPUTE RI, ZJ AND PR , QZ AT GRID POINT (IR , JZ)
!       ALSO, COMPUTE WEIGHTS WIJ FOR 4 CORNER GRID POINTS
!
         ri = rgrid(ir)
         zj = zgrid(jz)
         pr = (rad0 - ri)/delr
         qz = (zee0 - zj)/delz
         temp = pr*qz
         w22_bi(index1d,kin) = temp                   !p*q               
         w21_bi(index1d,kin) = pr - temp              !p*(1-q)
         w12_bi(index1d,kin) = qz - temp              !q*(1-p)
         w11_bi(index1d,kin) = 1 + temp - (pr + qz)   !(1-p)*(1-q)

      END DO
      END DO

      END SUBROUTINE setbivariate

      SUBROUTINE bivariate4pt (resp_rz, resp_su, kin)
      IMPLICIT NONE
      INTEGER, INTENT(in)      :: kin
      REAL(rprec), INTENT(in)  :: resp_rz(nrgrid*nzgrid)
      REAL(rprec), INTENT(out) :: resp_su(ndim)
      INTEGER :: jsu, irz11, irz12, irz21, irz22
!-----------------------------------------------
!
!     DETERMINE resp_uv ON THE S-U GRID BY 4 PT BIVARIATE INTERPOLATION
!     OF resp_rz WHICH IS ON AN R-Z MESH, ALL AT A FIXED TOROIDAL PLANE
!
!     USES 2-D INTERPOLATION BASED ON THE FOUR POINT FORMULA
!     IN ABRAMOWITZ AND STEGUN, EQ. 25.2.66 
!
!     resp_rz(1:nrgrid, 1:nzgrid) :   response on the rectangular R, Z grid
!
!     resp_su(1:ns,1:nu)          :   interpolated response on s,u grid (as 1d array)
!
!-----------------------------------------------
      DO jsu = 1, ndim
         
         irz11 = irz11_bi(jsu,kin)
         irz22 = irz22_bi(jsu,kin)
         irz12 = irz12_bi(jsu,kin)
         irz21 = irz21_bi(jsu,kin)
!
!       COMPUTE RESPONSE AT R, PHI (fixed kin), Z BY 4-PT INTERPOLATION
!
         resp_su(jsu) = w11_bi(jsu,kin)*resp_rz(irz11) 
     &                + w22_bi(jsu,kin)*resp_rz(irz22) 
     &                + w21_bi(jsu,kin)*resp_rz(irz21) 
     &                + w12_bi(jsu,kin)*resp_rz(irz12)

      END DO

      END SUBROUTINE bivariate4pt

      SUBROUTINE clear_bivar

      IF (ALLOCATED(irz11_bi)) DEALLOCATE (irz11_bi)
      IF (ALLOCATED(irz12_bi)) DEALLOCATE (irz12_bi)
      IF (ALLOCATED(irz21_bi)) DEALLOCATE (irz21_bi)
      IF (ALLOCATED(irz22_bi)) DEALLOCATE (irz22_bi)
      IF (ALLOCATED(w11_bi)) DEALLOCATE (w11_bi) 
      IF (ALLOCATED(w12_bi)) DEALLOCATE (w12_bi)
      IF (ALLOCATED(w21_bi)) DEALLOCATE (w21_bi)
      IF (ALLOCATED(w22_bi)) DEALLOCATE (w22_bi)
  
      END SUBROUTINE clear_bivar

      END MODULE bivariate
