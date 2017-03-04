

      MODULE read_boozer_mod
!
!     USAGE:
!
!     Use READ_BOOZ_MOD to include variables dynamically ALlocated
!     in the module
!     Call DEALLOCATE_READ_BOOZER to free this memory when it is no longer needed
!
      USE stel_kinds
      IMPLICIT NONE

C-----------------------------------------------
C   L O C A L   P A R A M E T E R S
C-----------------------------------------------
! Variable names (vn_...) : put eventually into library, used by read_wout too...
      CHARACTER(LEN=*), PARAMETER ::
     1 vn_nfp="nfp_b", vn_ns="ns_b", vn_aspect="aspect_b",
     2 vn_rmax="rmax_b", vn_rmin="rmin_b", vn_betaxis="betaxis_b",
     3 vn_mboz="mboz_b", vn_nboz="nboz_b", vn_mnboz="mnboz_b",
     4 vn_version="version", vn_iota="iota_b", vn_pres="pres_b",
     5 vn_beta="beta_b", vn_phip="phip_b", vn_phi="phi_b",
     6 vn_bvco="bvco_b", vn_buco="buco_b", vn_ixm="ixm_b",
     7 vn_ixn="ixn_b", vn_bmnc="bmnc_b", vn_rmnc="rmnc_b",
     8 vn_zmns="zmns_b", vn_pmns="pmns_b", vn_gmnc="gmn_b",
     9 vn_bmns="bmns_b", vn_rmns="rmns_b", vn_zmnc="zmnc_b", 
     A vn_pmnc="pmnc_b", vn_gmns="gmns_b", vn_lasym="lasym",
     B vn_jlist="jlist"


!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: mnboz_b, mboz_b, nboz_b, nfp_b, ns_b
      INTEGER, DIMENSION(:), ALLOCATABLE :: idx_b, ixm_b, ixn_b
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: iota_b, pres_b,
     1    phip_b, phi_b, beta_b, buco_b, bvco_b
      REAL(rprec) :: aspect_b, rmax_b, rmin_b, betaxis_b
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::
     1   bmnc_b, rmnc_b, zmns_b, pmns_b, gmnc_b, packed2d
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::
     1   bmns_b, rmns_b, zmnc_b, pmnc_b, gmns_b
      LOGICAL :: lasym_b = .FALSE.

      CONTAINS

      SUBROUTINE read_boozer_file (file_or_extension, ierr, iopen)
      USE safe_open_mod
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: ierr
      INTEGER, OPTIONAL :: iopen
      CHARACTER(LEN=*) :: file_or_extension
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: unit_booz = 14
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: iunit
      CHARACTER(len=LEN_TRIM(file_or_extension)+10) :: filename
      LOGICAL :: isnc
C-----------------------------------------------
!
!     THIS SUBROUTINE READS THE BOOZMN FILE CREATED BY THE BOOZ_XFORM CODE
!     AND STORES THE DATA IN THE READ_BOOZ_MOD MODULE
!
!     CHECK FOR netcdf FILE EXTENSION (*.nc)
!
      filename = 'boozmn'
      CALL parse_extension(filename, file_or_extension, isnc)

      IF (isnc) THEN

         CALL read_boozer_nc(filename, ierr)




      ELSE
         iunit = unit_booz
         CALL safe_open (iunit, ierr, filename, 'old', 'unformatted')
         IF (ierr .eq. 0) CALL read_boozer_bin(iunit, ierr)
         CLOSE(unit=iunit)
      END IF

      IF (PRESENT(iopen)) iopen = ierr

      END SUBROUTINE read_boozer_file


      SUBROUTINE read_boozer_nc(filename, ierr)
      USE stel_constants, ONLY: zero
      USE ezcdf
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: ierr
      CHARACTER(LEN=*) :: filename
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, DIMENSION(3) :: dimlens
      INTEGER :: nbooz, nsval, ilist
      INTEGER, ALLOCATABLE, DIMENSION(:) :: jlist
      CHARACTER(LEN=38) :: version
C-----------------------------------------------
! Open cdf File
      call cdf_open(nbooz,filename,'r', ierr)
      IF (ierr .ne. 0) THEN
         PRINT *,' Error opening boozmn .nc file'
         RETURN
      END IF

! Read in scalar variables
      CALL cdf_read(nbooz, vn_nfp, nfp_b)
      CALL cdf_read(nbooz, vn_ns, ns_b)
      CALL cdf_read(nbooz, vn_aspect, aspect_b)
      CALL cdf_read(nbooz, vn_rmax, rmax_b)
      CALL cdf_read(nbooz, vn_rmin, rmin_b)
      CALL cdf_read(nbooz, vn_betaxis, betaxis_b)
      CALL cdf_read(nbooz, vn_mboz, mboz_b)
      CALL cdf_read(nbooz, vn_nboz, nboz_b)
      CALL cdf_read(nbooz, vn_mnboz, mnboz_b)
      CALL cdf_read(nbooz, vn_version, version)
      CALL cdf_read(nbooz, vn_lasym, lasym_b)

!  1D arrays (skip inquiry statements for now, assume correct in file)
      IF (ALLOCATED(iota_b)) CALL read_boozer_deallocate
      ALLOCATE (iota_b(ns_b), pres_b(ns_b), beta_b(ns_b), phip_b(ns_b),
     1  phi_b(ns_b), bvco_b(ns_b), buco_b(ns_b), idx_b(ns_b), 
     2  ixm_b(mnboz_b), ixn_b(mnboz_b), stat=ierr)
      IF (ierr .ne. 0) THEN
        PRINT *,' Allocation error in read_boozer_file'
        RETURN
      END IF

      CALL cdf_read(nbooz, vn_iota, iota_b)
      CALL cdf_read(nbooz, vn_pres, pres_b)
      CALL cdf_read(nbooz, vn_beta, beta_b)
      CALL cdf_read(nbooz, vn_phip, phip_b)
      CALL cdf_read(nbooz, vn_phi,  phi_b)
      CALL cdf_read(nbooz, vn_bvco, bvco_b)
      CALL cdf_read(nbooz, vn_buco, buco_b)
      CALL cdf_read(nbooz, vn_ixm, ixm_b)
      CALL cdf_read(nbooz, vn_ixn, ixn_b)

      CALL cdf_inquire(nbooz, vn_jlist, dimlens)
      ALLOCATE (jlist(1:dimlens(1)), stat=ierr)
      CALL cdf_read(nbooz, vn_jlist, jlist)

      idx_b = 0
      ilist = SIZE(jlist)
      DO ilist = 1, SIZE(jlist)
         nsval = jlist(ilist)
         idx_b(nsval) = 1
      END DO

!  2D arrays
      ALLOCATE (bmnc_b(mnboz_b,ns_b), rmnc_b(mnboz_b,ns_b),
     1  zmns_b(mnboz_b,ns_b), pmns_b(mnboz_b,ns_b), 
     2  gmnc_b(mnboz_b,ns_b), packed2d(mnboz_b, ilist), stat=ierr)
      IF (ierr .ne. 0) THEN
         PRINT *,' Allocation error in read_boozer_file'
         RETURN
      END IF

!
!  Note: Must unpack these 2D arrays, only jlist-ed radial nodes store in file
!
      rmnc_b = 0; zmns_b = 0; pmns_b = 0; bmnc_b = 0; gmnc_b = 0
      CALL unpack_cdf(nbooz, vn_bmnc, bmnc_b)
      CALL unpack_cdf(nbooz, vn_rmnc, rmnc_b)
      CALL unpack_cdf(nbooz, vn_zmns, zmns_b)
      CALL unpack_cdf(nbooz, vn_pmns, pmns_b)
      CALL unpack_cdf(nbooz, vn_gmnc, gmnc_b)

      IF (lasym_b) THEN
      ALLOCATE (bmns_b(mnboz_b,ns_b), rmns_b(mnboz_b,ns_b),
     1  zmnc_b(mnboz_b,ns_b), pmnc_b(mnboz_b,ns_b),
     2  gmns_b(mnboz_b,ns_b), stat=ierr)
      IF (ierr .ne. 0) THEN
         PRINT *,' Allocation error in read_boozer_file'
         RETURN
      END IF
      rmns_b = 0; zmnc_b = 0; pmnc_b = 0; bmns_b = 0; gmns_b = 0
      CALL unpack_cdf(nbooz, vn_bmns, bmns_b)
      CALL unpack_cdf(nbooz, vn_rmns, rmns_b)
      CALL unpack_cdf(nbooz, vn_zmnc, zmnc_b)
      CALL unpack_cdf(nbooz, vn_pmnc, pmnc_b)
      CALL unpack_cdf(nbooz, vn_gmns, gmns_b)
      END IF


      DEALLOCATE (jlist, packed2d)

! Close cdf File
      CALL cdf_close(nbooz, ierr)

      END SUBROUTINE read_boozer_nc


      SUBROUTINE read_boozer_bin(iunit, ierr)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: ierr, iunit
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: nsval, jsize, js
      CHARACTER(LEN=38) :: version
C-----------------------------------------------

      READ(iunit, iostat=ierr, err=100) nfp_b, ns_b, aspect_b,
     1   rmax_b, rmin_b, betaxis_b

      IF (ALLOCATED(iota_b)) CALL read_boozer_deallocate
      ALLOCATE (iota_b(ns_b), pres_b(ns_b), beta_b(ns_b), phip_b(ns_b),
     1  phi_b(ns_b), bvco_b(ns_b), buco_b(ns_b), idx_b(ns_b), stat=ierr)
      IF (ierr .ne. 0) THEN
        PRINT *,' Allocation error in read_boozer_file'
        RETURN
      END IF
 
      iota_b(1) = 0; pres_b(1) = 0; beta_b(1) = 0
      phip_b(1) = 0; phi_b(1) = 0; bvco_b(1) = 0
      buco_b(1) = 0

      DO nsval = 2, ns_b
         READ(iunit, iostat=ierr, err=100) iota_b(nsval),
     1   pres_b(nsval), beta_b(nsval), phip_b(nsval), phi_b(nsval),
     2   bvco_b(nsval), buco_b(nsval)
      END DO

      READ(iunit, iostat=ierr, err=100) mboz_b, nboz_b, mnboz_b, jsize
      READ(iunit, iostat=js) version, lasym_b

      ALLOCATE (bmnc_b(mnboz_b,ns_b), rmnc_b(mnboz_b,ns_b),
     1  zmns_b(mnboz_b,ns_b), pmns_b(mnboz_b,ns_b), 
     2  gmnc_b(mnboz_b,ns_b), ixm_b(mnboz_b), ixn_b(mnboz_b), stat=ierr)
      IF (ierr .ne. 0) THEN
         PRINT *,' Allocation error in read_boozer_file'
         RETURN
      END IF

      idx_b = 0
      ixm_b = 0
      rmnc_b = 0; zmns_b = 0; pmns_b = 0; bmnc_b = 0; gmnc_b = 0

      IF (lasym_b) THEN
      ALLOCATE (bmns_b(mnboz_b,ns_b), rmns_b(mnboz_b,ns_b),
     1  zmnc_b(mnboz_b,ns_b), pmnc_b(mnboz_b,ns_b), 
     2  gmns_b(mnboz_b,ns_b), stat=ierr)
      IF (ierr .ne. 0) THEN
         PRINT *,' Allocation error in read_boozer_file'
         RETURN
      END IF

      rmns_b = 0; zmnc_b = 0; pmnc_b = 0; bmns_b = 0; gmns_b = 0

      END IF

!
!     idx_b:  = 0, surface data not in file; = 1, surface data in file
!
      READ(iunit,iostat=ierr,err=100) ixn_b(:mnboz_b), ixm_b(:mnboz_b)

      DO js = 1, jsize
        READ(iunit, iostat=ierr, END=200, err=100) nsval
        IF ((nsval.gt.ns_b) .or. (nsval.le.0)) CYCLE

        idx_b(nsval) = 1

        READ(iunit,iostat=ierr,err=100, END=200) bmnc_b(:mnboz_b,nsval),
     1       rmnc_b(:mnboz_b,nsval), zmns_b(:mnboz_b,nsval), 
     2       pmns_b(:mnboz_b,nsval), gmnc_b(:mnboz_b,nsval)

        IF (.not.lasym_b) CYCLE

        READ(iunit,iostat=ierr,err=100, END=200) bmns_b(:mnboz_b,nsval),
     1       rmns_b(:mnboz_b,nsval), zmnc_b(:mnboz_b,nsval), 
     2       pmnc_b(:mnboz_b,nsval), gmns_b(:mnboz_b,nsval)
      END DO

 100  CONTINUE
      IF (ierr .gt. 0)
     1    PRINT *,' Error reading in subroutine read_boozer_file:',
     2            ' ierr = ', ierr
 200  CONTINUE
      IF (ierr .lt. 0) ierr = 0       !End-of-file, ok
      CLOSE(iunit)

      END SUBROUTINE read_boozer_bin


      SUBROUTINE read_boozer_deallocate
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: istat
!-----------------------------------------------

      IF (ALLOCATED(iota_b)) DEALLOCATE (iota_b, pres_b, beta_b,
     1    phip_b, phi_b, bvco_b, buco_b, idx_b, stat = istat)

      IF (ALLOCATED(bmnc_b)) DEALLOCATE (bmnc_b, rmnc_b,
     1   zmns_b, pmns_b, gmnc_b, ixm_b, ixn_b, stat = istat)

      IF (ALLOCATED(bmns_b)) DEALLOCATE (bmns_b, rmns_b,
     1   zmnc_b, pmnc_b, gmns_b, stat = istat)

      END SUBROUTINE read_boozer_deallocate



      SUBROUTINE unpack_cdf (nbooz, var_name, array2d)
      USE stel_kinds, ONLY: rprec
      USE ezcdf
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nbooz
      INTEGER :: nsval, icount
      REAL(rprec), DIMENSION(:,:), INTENT(out) :: array2d
      CHARACTER(LEN=*), INTENT(in) :: var_name
C-----------------------------------------------
!
!     Read into temporary packed array, packed2d
!
      CALL cdf_read(nbooz, var_name, packed2d)

      array2d = 0; icount = 1

      DO nsval = 1, ns_b
         IF (idx_b(nsval) .eq. 1) THEN
            array2d(:,nsval) = packed2d(:,icount)
            icount = icount + 1
         END IF
      END DO

      END SUBROUTINE unpack_cdf



      SUBROUTINE write_boozer_nc(extension, ierr)
      USE stel_constants, ONLY: zero
!DEC$ IF DEFINED (NETCDF)
      USE ezcdf
!DEC$ ENDIF
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: ierr
      CHARACTER(LEN=*) :: extension
!DEC$ IF DEFINED (NETCDF)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, DIMENSION(3) :: dimlens
      INTEGER :: nbooz, nsval, ilist, jsize, i
      INTEGER, ALLOCATABLE, DIMENSION(:) :: jlist
      CHARACTER(LEN=*), PARAMETER, DIMENSION(1) ::
     1             r1dim = (/'radius'/), mn1dim = (/'mn_mode'/),
     2             j1dim = (/'comput_surfs'/)
      CHARACTER(LEN=*), DIMENSION(2), PARAMETER :: 
     1             r2dim = (/'mn_modes','pack_rad'/)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      CHARACTER(LEN=*), PARAMETER :: version =
     1   "Boozer Transformation Code Version 2.0"
C-----------------------------------------------
! Open cdf File
      CALL cdf_open(nbooz,'boozmn_' // TRIM(extension) // '.nc',
     1             'w', ierr)
      IF (ierr .ne. 0) THEN
         PRINT *,' Error opening boozmn_'// TRIM(extension) 
     1           // '.nc file'
         RETURN
      END IF

      jsize = COUNT(idx_b == 1)
      i = 1
      ALLOCATE(jlist(jsize),packed2d(mnboz_b,jsize))
      DO ilist = 1, ns_b
         IF (idx_b(ilist) .le. 0) CYCLE
         jlist(i) = ilist
         i = i + 1
      END DO
! Define Variables
! Scalars
      CALL cdf_define(nbooz, vn_nfp, nfp_b)
      CALL cdf_define(nbooz, vn_ns, ns_b)
      CALL cdf_define(nbooz, vn_aspect, aspect_b)
      CALL cdf_define(nbooz, vn_rmax, rmax_b)
      CALL cdf_define(nbooz, vn_rmin, rmin_b)
      CALL cdf_define(nbooz, vn_betaxis, betaxis_b)
      CALL cdf_define(nbooz, vn_mboz, mboz_b)
      CALL cdf_define(nbooz, vn_nboz, nboz_b)
      CALL cdf_define(nbooz, vn_mnboz, mnboz_b)
      CALL cdf_define(nbooz, vn_version, version)
      CALL cdf_define(nbooz, vn_lasym, lasym_b)
! 1D Arrays
      CALL cdf_define(nbooz, vn_iota, iota_b, dimname=r1dim)
      CALL cdf_define(nbooz, vn_pres, pres_b, dimname=r1dim)
      CALL cdf_define(nbooz, vn_beta, beta_b, dimname=r1dim)
      CALL cdf_define(nbooz, vn_phip, phip_b, dimname=r1dim)
      CALL cdf_define(nbooz, vn_phi,  phi_b, dimname=r1dim)
      CALL cdf_define(nbooz, vn_bvco, bvco_b, dimname=r1dim)
      CALL cdf_define(nbooz, vn_buco, buco_b, dimname=r1dim)
      CALL cdf_define(nbooz, vn_jlist, jlist, dimname=j1dim)
      CALL cdf_define(nbooz, vn_ixm, ixm_b, dimname=mn1dim)
      CALL cdf_define(nbooz, vn_ixn, ixn_b, dimname=mn1dim)
! 2D Arrays
      CALL cdf_define(nbooz, vn_bmnc, packed2d, dimname=r2dim)
      CALL cdf_define(nbooz, vn_rmnc, packed2d, dimname=r2dim)
      CALL cdf_define(nbooz, vn_zmns, packed2d, dimname=r2dim)
      CALL cdf_define(nbooz, vn_pmns, packed2d, dimname=r2dim)
      CALL cdf_define(nbooz, vn_gmnc, packed2d, dimname=r2dim)
      IF (lasym_b) THEN
      CALL cdf_define(nbooz, vn_bmns, packed2d, dimname=r2dim)
      CALL cdf_define(nbooz, vn_rmns, packed2d, dimname=r2dim)
      CALL cdf_define(nbooz, vn_zmnc, packed2d, dimname=r2dim)
      CALL cdf_define(nbooz, vn_pmnc, packed2d, dimname=r2dim)
      CALL cdf_define(nbooz, vn_gmns, packed2d, dimname=r2dim)
      ENDIF

! Write out scalars
      CALL cdf_write(nbooz, vn_nfp, nfp_b)
      CALL cdf_write(nbooz, vn_ns, ns_b)
      CALL cdf_write(nbooz, vn_aspect, aspect_b)
      CALL cdf_write(nbooz, vn_rmax, rmax_b)
      CALL cdf_write(nbooz, vn_rmin, rmin_b)
      CALL cdf_write(nbooz, vn_betaxis, betaxis_b)
      CALL cdf_write(nbooz, vn_mboz, mboz_b)
      CALL cdf_write(nbooz, vn_nboz, nboz_b)
      CALL cdf_write(nbooz, vn_mnboz, mnboz_b)
      CALL cdf_write(nbooz, vn_version, version)

!  1D arrays 

      CALL cdf_write(nbooz, vn_iota, iota_b)
      CALL cdf_write(nbooz, vn_pres, pres_b)
      CALL cdf_write(nbooz, vn_beta, beta_b)
      CALL cdf_write(nbooz, vn_phip, phip_b)
      CALL cdf_write(nbooz, vn_phi,  phi_b)
      CALL cdf_write(nbooz, vn_bvco, bvco_b)
      CALL cdf_write(nbooz, vn_buco, buco_b)
      CALL cdf_write(nbooz, vn_jlist, jlist)
      CALL cdf_write(nbooz, vn_ixm, ixm_b)
      CALL cdf_write(nbooz, vn_ixn, ixn_b)

!  Write packed 2D arrays
      CALL pack_cdf(nbooz, vn_bmnc,bmnc_b)
      CALL pack_cdf(nbooz, vn_rmnc,rmnc_b)
      CALL pack_cdf(nbooz, vn_zmns,zmns_b)
      CALL pack_cdf(nbooz, vn_zmns,pmns_b)
      CALL pack_cdf(nbooz, vn_rmnc,gmnc_b)
      !CALL cdf_write(nbooz, vn_bmnc, bmnc_b(:,2:ns_b))
      !CALL cdf_write(nbooz, vn_rmnc, rmnc_b(:,2:ns_b))
      !CALL cdf_write(nbooz, vn_zmns, zmns_b(:,2:ns_b))
      !CALL cdf_write(nbooz, vn_pmns, pmns_b(:,2:ns_b))
      !CALL cdf_write(nbooz, vn_gmnc, gmnc_b(:,2:ns_b))
      CALL cdf_write(nbooz, vn_lasym, lasym_b)
      IF (lasym_b) THEN
         CALL pack_cdf(nbooz, vn_bmns,bmnc_b)
         CALL pack_cdf(nbooz, vn_rmns,rmnc_b)
         CALL pack_cdf(nbooz, vn_zmnc,zmns_b)
         CALL pack_cdf(nbooz, vn_zmnc,pmns_b)
         CALL pack_cdf(nbooz, vn_rmns,gmnc_b)
         !CALL cdf_write(nbooz, vn_bmns, bmns_b(:,2:ns_b))
         !CALL cdf_write(nbooz, vn_rmns, rmns_b(:,2:ns_b))
         !CALL cdf_write(nbooz, vn_zmnc, zmnc_b(:,2:ns_b))
         !CALL cdf_write(nbooz, vn_pmnc, pmnc_b(:,2:ns_b))
         !CALL cdf_write(nbooz, vn_gmns, gmns_b(:,2:ns_b))
      ENDIF

      DEALLOCATE(jlist,packed2d)
! Close cdf File
      CALL cdf_close(nbooz, ierr)
!DEC$ ENDIF
      RETURN
      END SUBROUTINE write_boozer_nc
      


      SUBROUTINE write_boozer_bin(extension, ierr)
      USE stel_constants, ONLY: zero
      USE safe_open_mod
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: ierr
      CHARACTER(LEN=*) :: extension
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, iunit, istat, js, jsize
      REAL(rprec) :: version = 2.0

      CALL safe_open (iunit, istat, 'boozmn.' // extension, 'replace',
     1     'unformatted')
      IF (istat .ne. 0) THEN
         PRINT *,' istat = ', istat
         STOP 'Error opening boozmn file in XBOOZ_XFORM!'
      END IF

!
!     Write out surface quantities needed by ballooning code
!
      WRITE(iunit, iostat=istat, err=100)
     1   nfp_b, ns_b, aspect_b, rmax_b, rmin_b, betaxis_b

      DO js = 2, ns_b
         WRITE(iunit, iostat=istat, err=100) iota_b(js), pres_b(js),
     1   beta_b(js), phip_b(js), phi_b(js), bvco_b(js), buco_b(js)
      END DO

      jsize = COUNT(idx_b == 1)
!     SPH 070909: ADDED lasym to dump
      WRITE (iunit, iostat=istat, err=100) mboz_b, nboz_b, mnboz_b, 
     1                                     jsize
      WRITE (iunit, iostat=istat, err=100) version, lasym_b

      WRITE (iunit, iostat=istat, err=100) ixn_b(:mnboz_b), 
     1                                     ixm_b(:mnboz_b)
!
!     Write packed (in radial coordinate) 2D arrays
!
      DO i = 1, ns_b
        js = idx_b(i)
        IF (js.le.0 .or. js.gt.ns_b) CYCLE
        WRITE (iunit, iostat=istat, err=100) js
        WRITE (iunit, iostat=istat, err=100) bmnc_b(:mnboz_b,i),
     1       rmnc_b(:mnboz_b,i), zmns_b(:mnboz_b,i), 
     2       pmns_b(:mnboz_b,i), gmnc_b(:mnboz_b,i)
!SPH070909: WRITE OUT ASYMMETRIC PARTS
        IF (lasym_b) THEN
        WRITE (iunit, iostat=istat, err=100) bmnc_b(:mnboz_b,i),
     1       rmns_b(:mnboz_b,i), zmnc_b(:mnboz_b,i),
     2       pmnc_b(:mnboz_b,i), gmns_b(:mnboz_b,i)
        ENDIF	
      END DO


 100  CONTINUE
      IF (istat .gt. 0)
     1    PRINT *,' Error writing in subroutine write_boozmn:',
     2            ' istat = ', istat

      CLOSE(iunit)

      END SUBROUTINE write_boozer_bin



      SUBROUTINE pack_cdf (nbooz, var_name, array2d)
      USE stel_kinds, ONLY: rprec
      USE ezcdf
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nbooz
      INTEGER :: nsval, icount
      REAL(rprec), DIMENSION(:,:), INTENT(out) :: array2d
      CHARACTER(LEN=*), INTENT(in) :: var_name
C-----------------------------------------------
!
!     Read into temporary packed array, packed2d
!

      packed2d = 0; icount = 1

      DO nsval = 1, ns_b
         IF (idx_b(nsval) .eq. 1) THEN
            packed2d(:,icount) = array2d(:,nsval)
            icount = icount + 1
         END IF
      END DO
      CALL cdf_write(nbooz, var_name, packed2d)

      END SUBROUTINE pack_cdf


      END MODULE read_boozer_mod


 
