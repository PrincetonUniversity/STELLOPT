      MODULE write_mgrid
      USE biotsavart
      USE mgrid_mod, ONLY: vn_nextcur, vn_mgmode, vn_ir,
     1    vn_jz, vn_kp, vn_nfp, vn_rmin, vn_rmax,
     2    vn_zmin, vn_zmax, vn_coilgrp, vn_coilcur, 
     3    vn_br0, vn_bz0, vn_bp0
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: one = 1
      CHARACTER(LEN=*), PARAMETER ::
     1   coildim(2) = (/'stringsize          ',
     2                  'external_coil_groups'/),
     3   groupdim(1)= (/'external_coils'/),
     4   cylcoord(3)= (/'rad','zee','phi'/)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER     :: ir = 121, jz = 121
      INTEGER     :: iextc, nfp, nextcur, istat
      INTEGER     :: kp, kp2, jz2, kp_odd, jz_odd
      REAL(rprec), ALLOCATABLE, DIMENSION(:, :, :) :: br, bz, bp
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: extcur
      REAL(rprec) :: rmin, rmax, zmin, zmax
      REAL(rprec) :: fperiod, delr, delz, delp
      LOGICAL     :: lstell_sym=.false.
      CHARACTER(LEN=1) :: mgrid_mode='S'
      CHARACTER(LEN=70) :: mgrid_file, coil_file
      CHARACTER(LEN=60) :: mgrid_ext

      PRIVATE :: istat
C-----------------------------------------------

      CONTAINS
#if defined(NETCDF)
      SUBROUTINE write_mgrid_nc
      USE ezcdf
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ngrid, ig
      CHARACTER(LEN=100) :: temp
      CHARACTER(LEN=100), ALLOCATABLE, DIMENSION(:) :: 
     1                    vn_br, vn_bp, vn_bz
C-----------------------------------------------
      mgrid_file = TRIM(mgrid_file) // '.nc'

      CALL cdf_open(ngrid,mgrid_file,'w',istat)
      if (istat .ne. 0) STOP 'Error opening WOUT.nc file VMEC WROUT'

!
!     DEFINE DATA VARIABLES, DIMENSION NAMES
!
      CALL cdf_define(ngrid, vn_ir, ir)
      CALL cdf_define(ngrid, vn_jz, jz)
      CALL cdf_define(ngrid, vn_kp, kp)
      CALL cdf_define(ngrid, vn_nfp, nfp)
      CALL cdf_define(ngrid, vn_nextcur, nextcur)
      CALL cdf_define(ngrid, vn_rmin, rmin)
      CALL cdf_define(ngrid, vn_zmin, zmin)
      CALL cdf_define(ngrid, vn_rmax, rmax)
      CALL cdf_define(ngrid, vn_zmax, zmax)
      IF (nextcur .eq. 1) THEN
         CALL cdf_define(ngrid, vn_coilgrp,coil_group(1)%s_name) 
      ELSE
         CALL cdf_define(ngrid, vn_coilgrp,coil_group(1:nextcur)%s_name,
     1                   dimname=coildim)
      END IF
      CALL cdf_define(ngrid, vn_mgmode, mgrid_mode)
      CALL cdf_define(ngrid, vn_coilcur, extcur(1:nextcur),
     1                dimname=groupdim)
!
!     STORED AS 3D ARRAYS (ACTUALLY 4D, BUT CUT THROUGH IG)
!
      ALLOCATE (vn_br(nextcur), vn_bz(nextcur), vn_bp(nextcur))

      DO ig = 1, nextcur
         write (temp, '(a,i3.3)') "_",ig
         vn_br(ig) = vn_br0 // temp
         vn_bp(ig) = vn_bp0 // temp
         vn_bz(ig) = vn_bz0 // temp
         CALL cdf_define(ngrid, vn_br(ig), br, dimname=cylcoord)
         CALL cdf_define(ngrid, vn_bp(ig), bp, dimname=cylcoord)
         CALL cdf_define(ngrid, vn_bz(ig), bz, dimname=cylcoord)
      END DO


!
!     WRITE OUT DATA
!
      CALL cdf_write(ngrid, vn_ir, ir)
      CALL cdf_write(ngrid, vn_jz, jz)
      CALL cdf_write(ngrid, vn_kp, kp)
      CALL cdf_write(ngrid, vn_nfp, nfp)
      CALL cdf_write(ngrid, vn_nextcur, nextcur)
      CALL cdf_write(ngrid, vn_rmin, rmin)
      CALL cdf_write(ngrid, vn_zmin, zmin)
      CALL cdf_write(ngrid, vn_rmax, rmax)
      CALL cdf_write(ngrid, vn_zmax, zmax)
      IF (nextcur .eq. 1) THEN
         CALL cdf_write(ngrid, vn_coilgrp, coil_group(1)%s_name)
      ELSE
         CALL cdf_write(ngrid, vn_coilgrp, coil_group(1:nextcur)%s_name)
      END IF

!
!     SET UP CYLINDRICAL COMPONENTS OF MAGNETIC FIELD ON GRID
!     SUM OVER CURRENT GROUPS IG = 1,NEXTCUR
!     NOTE: USER MUST SUPPLY SUBROUTINE "BFIELD" TO COMPUTE THESE VALUES
!
      GROUPS: DO ig = 1,nextcur

         CALL compute_bfield(ig)

         CALL cdf_write(ngrid, vn_br(ig), br)
         CALL cdf_write(ngrid, vn_bp(ig), bp)
         CALL cdf_write(ngrid, vn_bz(ig), bz)

      END DO GROUPS 

      CALL cdf_write(ngrid, vn_mgmode, mgrid_mode)
      CALL cdf_write(ngrid, vn_coilcur, extcur(1:nextcur))

      CALL cdf_close(ngrid)

      DEALLOCATE (vn_br, vn_bz, vn_bp)

      END SUBROUTINE write_mgrid_nc
#endif

      SUBROUTINE write_mgrid_bin
      USE safe_open_mod
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: igrid0 = 10
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER  :: igrid, ig, i
C-----------------------------------------------

      igrid = igrid0
      mgrid_file = TRIM(mgrid_file) // '.bin'
      CALL safe_open (igrid, istat, TRIM(mgrid_file),
     1   'replace', 'unformatted')
      
      IF (istat .ne. 0) THEN
         WRITE (6,*) 'XGRID could not create ', TRIM(mgrid_file)
         WRITE (6,*) 'IOSTAT = ', istat,' IUNIT = ', igrid
         STOP 20
      END IF

!
!     MARK NEW-STYLE CODE WITH NEXTCUR < 0 SO READ_MGRID CAN DISTINGUISH
!     FROM OLD-STYLE MGRID FILE
!
      WRITE(igrid) ir, jz, kp, nfp, -nextcur
      WRITE(igrid) rmin, zmin, rmax, zmax
      WRITE(igrid) (coil_group(i) % s_name, i=1,nextcur)

!
!     SET UP CYLINDRICAL COMPONENTS OF MAGNETIC FIELD ON GRID
!     SUM OVER CURRENT GROUPS IG = 1,NEXTCUR
!     NOTE: USER MUST SUPPLY SUBROUTINE "BFIELD" TO COMPUTE THESE VALUES
!
      GROUPS: DO ig = 1,nextcur

         CALL compute_bfield (ig)


         WRITE (igrid, iostat=istat) br, bp, bz
!        WRITE(igrid,'(1p,9e12.4)') (((br(i,j,k), bz(i,j,k),              !OLD STYLE
!    1         bp(i,j,k), i=1,ir), j=1,jz),k=1,kp)
         IF (istat .ne. 0) STOP ' Error writing bfield components'

      END DO GROUPS 

      WRITE (igrid) mgrid_mode
      WRITE (igrid) extcur(1:nextcur)

!
!     ADD RECONSTRUCTION/OBSERVATION STUFF HERE
!

      CLOSE (igrid)

      END SUBROUTINE write_mgrid_bin

       
      SUBROUTINE compute_bfield (ig)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ig
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: icoil, numcoils, k, j, i, numfils, curindex(1)
      REAL(rprec) :: rad, zee, phi, extcur_ig, extcur_ig_inv
C-----------------------------------------------
      numcoils = coil_group(ig) % ncoil
      curindex = MAXLOC(abs(coil_group(ig)%coils(1:numcoils)%current))
      extcur_ig = coil_group(ig)%coils(curindex(1))%current
      
      extcur(ig) = extcur_ig
      IF (extcur_ig .eq. zero) WRITE (6, 100) ig
 100  FORMAT ('In COILS file, current(s) vanished for coil group ',i4)

      numfils = 0
      DO icoil = 1, numcoils
         numfils = numfils +                                                   &
     &      MAX(SIZE(coil_group(ig) % coils(icoil) % xnod,2),2) - 1
!
!        Compute field for unit current
!
!  JDH 2011-08-16. Comment out below - don't want to change bsc_coil data
!    Adjust B values late on with Raw test
!    Will behave differently if extcur_ig == 0.
!         IF (extcur_ig .ne. zero) THEN
!            coil_group(ig) % coils(icoil) % current =
!     1      coil_group(ig) % coils(icoil) % current /extcur_ig
!         ELSE
!            coil_group(ig) % coils(icoil) % current = 1
!         END IF
      END DO

      WRITE (6, '(/,2a)')          ' COIL GROUP          :  ', 
     1      TRIM(coil_group(ig) % s_name)
      WRITE (6, '(a,i6,a,i6)') ' TOTAL COILS IN GROUP: ', numcoils,
     1' TOTAL FILAMENTS: ', numfils

      k = 1                     ! this is always a symmetry plane
      phi = (k-1)*delp

!$omp parallel do private(zee, rad)
      DO j=1,jz2 + jz_odd
         zee = zmin + (j-1)*delz
         DO i=1,ir
            rad = rmin + (i-1)*delr
            CALL bfield (rad, phi, zee, br(i,j,k),
     1                   bp(i,j,k), bz(i,j,k), ig)

            IF (lstell_sym) THEN
               br(i,jz+1-j,k) = -br(i,j,k)
               bz(i,jz+1-j,k) =  bz(i,j,k)
               bp(i,jz+1-j,k) =  bp(i,j,k)
            END IF
         END DO
      END DO
!$omp end parallel do
      IF (kp .ne. 1) 
     1   WRITE (6, '(a,i4,a,i4,a)') ' K = ',k,' (OUT OF ',KP,')'

!$omp parallel do private(zee, rad, phi)
      DO k=2,kp2+kp_odd
         phi = (k-1)*delp
         DO j=1,jz
            zee = zmin + (j-1)*delz
            DO i=1,ir
               rad = rmin + (i-1)*delr
               CALL bfield (rad, phi, zee, br(i,j,k),
     1                      bp(i,j,k), bz(i,j,k), ig)
               IF (lstell_sym) THEN
                  br(i,jz+1-j,kp+2-k) = -br(i,j,k)
                  bz(i,jz+1-j,kp+2-k) =  bz(i,j,k)
                  bp(i,jz+1-j,kp+2-k) =  bp(i,j,k)
               END IF
            END DO
         END DO
         WRITE (6,'(a,i4)') ' K = ',k
      END DO
!$omp end parallel do

      IF ((kp_odd == 0) .and. lstell_sym) THEN       ! another symmetry plane
         k = kp2 + 1
         phi = (k-1)*delp
!$omp parallel do private(zee, rad)
         DO j=1,jz2 + jz_odd
            zee = zmin + (j-1)*delz
            DO i=1,ir
               rad = rmin + (i-1)*delr
               CALL bfield (rad, phi, zee, br(i,j,k),
     1                      bp(i,j,k), bz(i,j,k), ig)
               br(i,jz+1-j,k) = -br(i,j,k)
               bz(i,jz+1-j,k) =  bz(i,j,k)
               bp(i,jz+1-j,k) =  bp(i,j,k)
            END DO
         END DO
!$omp end parallel do
         WRITE (6,'(a,i4)') ' K = ',k
      END IF

      IF (mgrid_mode .eq. 'R') THEN
!  JDH 2011-080-16. Comment out below (1 line)
!         br = br*extcur_ig;  bp = bp*extcur_ig; bz = bz*extcur_ig    ! "Raw" fields
         IF (ig .lt. 10) WRITE (iextc, 220) ig, extcur_ig, extcur_ig      ! currents for "raw" mode
         IF (ig .ge. 10) WRITE (iextc, 225) ig, extcur_ig, extcur_ig
      ELSE
!  JDH 2011-08-16. Add below, scale b. 2011-08-24 - clean up extcur_ig_inv 
         IF (ABS(extcur_ig) .gt. 1.E-100_rprec) THEN
            extcur_ig_inv = 1. / extcur_ig
         ELSE
            extcur_ig_inv = 1.E100_rprec
         ENDIF
!         extcur_ig_inv = 1. / MAX(extcur_ig, 1.E-100)
         br = br*extcur_ig_inv;  bp = bp*extcur_ig_inv; 
         bz = bz*extcur_ig_inv                                 ! "Scaled" fields
!  JDH 2011-08-16. end addition 
         IF (ig .lt. 10) WRITE (iextc, 210) ig, extcur_ig      ! currents for "scaled" mode
         IF (ig .ge. 10) WRITE (iextc, 215) ig, extcur_ig
      END IF


 210  FORMAT('EXTCUR(', i1,')  = ', 1p,e22.14)
 215  FORMAT('EXTCUR(', i2,') = ', 1p,e22.14)
 220  FORMAT('EXTCUR(', i1,')  = ', 1p,e22.14,'/',e22.14)
 225  FORMAT('EXTCUR(', i2,') = ', 1p,e22.14,'/',e22.14)

      END SUBROUTINE compute_bfield

      END MODULE write_mgrid
