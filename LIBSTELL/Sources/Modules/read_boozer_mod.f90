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
#if defined(NETCDF)
!-----------------------------------------------
!   L O C A L   P A R A M E T E R S
!-----------------------------------------------
! Variable names (vn_...) : put eventually into library, used by read_wout too...
      CHARACTER(LEN=*), PARAMETER ::                                    &
       vn_nfp="nfp_b", vn_ns="ns_b", vn_aspect="aspect_b",              &
       vn_rmax="rmax_b", vn_rmin="rmin_b", vn_betaxis="betaxis_b",      &
       vn_mboz="mboz_b", vn_nboz="nboz_b", vn_mnboz="mnboz_b",          &
       vn_version="version", vn_iota="iota_b", vn_pres="pres_b",        &
       vn_beta="beta_b", vn_phip="phip_b", vn_phi="phi_b",              &
       vn_bvco="bvco_b", vn_buco="buco_b", vn_ixm="ixm_b",              &
       vn_ixn="ixn_b", vn_bmnc="bmnc_b", vn_rmnc="rmnc_b",              &
       vn_zmns="zmns_b", vn_pmns="pmns_b", vn_gmnc="gmn_b",             &
       vn_bmns="bmns_b", vn_rmns="rmns_b", vn_zmnc="zmnc_b",            &
       vn_pmnc="pmnc_b", vn_gmns="gmns_b", vn_lasym="lasym",            &
       vn_jlist="jlist"
#endif
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: mnboz_b, mboz_b, nboz_b, nfp_b, ns_b
      INTEGER, DIMENSION(:), ALLOCATABLE :: idx_b, ixm_b, ixn_b
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: iota_b, pres_b,         &
         phip_b, phi_b, beta_b, buco_b, bvco_b
      REAL(rprec) :: aspect_b, rmax_b, rmin_b, betaxis_b
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::                       &
        bmnc_b, rmnc_b, zmns_b, pmns_b, gmnc_b, packed2d
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::                       &
         bmns_b, rmns_b, zmnc_b, pmnc_b, gmns_b
      LOGICAL :: lasym_b = .FALSE.

      CONTAINS

      SUBROUTINE read_boozer_file (file_or_extension, ierr, iopen)
      USE safe_open_mod
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: ierr
      INTEGER, OPTIONAL :: iopen
      CHARACTER(LEN=*) :: file_or_extension
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: unit_booz = 14
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: iunit
      CHARACTER(len=LEN_TRIM(file_or_extension)+10) :: filename
      LOGICAL :: isnc
!-----------------------------------------------
!
!     THIS SUBROUTINE READS THE BOOZMN FILE CREATED BY THE BOOZ_XFORM CODE
!     AND STORES THE DATA IN THE READ_BOOZ_MOD MODULE
!
!     CHECK FOR netcdf FILE EXTENSION (*.nc)
!
      filename = 'boozmn'
      CALL parse_extension(filename, file_or_extension, isnc)

      IF (isnc) THEN
#if defined(NETCDF)
         CALL read_boozer_nc(filename, ierr)
#else
         PRINT *, "NETCDF wout file can not be opened on this platform"
         ierr = -100
#endif
      ELSE
         iunit = unit_booz
         CALL safe_open (iunit, ierr, filename, 'old', 'unformatted')
         IF (ierr .eq. 0) CALL read_boozer_bin(iunit, ierr)
         CLOSE(unit=iunit)
      END IF

      IF (PRESENT(iopen)) iopen = ierr

      END SUBROUTINE read_boozer_file

#if defined(NETCDF)
      SUBROUTINE read_boozer_nc(filename, ierr)
      USE stel_constants, ONLY: zero
      USE ezcdf
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: ierr
      CHARACTER(LEN=*) :: filename
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(3) :: dimlens
      INTEGER :: nbooz, nsval, ilist
      INTEGER, ALLOCATABLE, DIMENSION(:) :: jlist
      CHARACTER(LEN=38) :: version
!-----------------------------------------------
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
      ALLOCATE (iota_b(ns_b), pres_b(ns_b), beta_b(ns_b), phip_b(ns_b), &
        phi_b(ns_b), bvco_b(ns_b), buco_b(ns_b), idx_b(ns_b),           &
        ixm_b(mnboz_b), ixn_b(mnboz_b), stat=ierr)
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
      ALLOCATE (bmnc_b(mnboz_b,ns_b), rmnc_b(mnboz_b,ns_b),             &
        zmns_b(mnboz_b,ns_b), pmns_b(mnboz_b,ns_b),                     &
        gmnc_b(mnboz_b,ns_b), packed2d(mnboz_b, ilist), stat=ierr)
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
      ALLOCATE (bmns_b(mnboz_b,ns_b), rmns_b(mnboz_b,ns_b),             &
        zmnc_b(mnboz_b,ns_b), pmnc_b(mnboz_b,ns_b),                     &
        gmns_b(mnboz_b,ns_b), stat=ierr)
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
#endif

      SUBROUTINE read_boozer_bin(iunit, ierr)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: ierr, iunit
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: nsval, jsize, js
      CHARACTER(LEN=38) :: version
!-----------------------------------------------

      READ(iunit, iostat=ierr, err=100) nfp_b, ns_b, aspect_b,          &
           rmax_b, rmin_b, betaxis_b

      IF (ALLOCATED(iota_b)) CALL read_boozer_deallocate
      ALLOCATE (iota_b(ns_b), pres_b(ns_b), beta_b(ns_b), phip_b(ns_b),   &
        phi_b(ns_b), bvco_b(ns_b), buco_b(ns_b), idx_b(ns_b), stat=ierr)
      IF (ierr .ne. 0) THEN
        PRINT *,' Allocation error in read_boozer_file'
        RETURN
      END IF
 
      iota_b(1) = 0; pres_b(1) = 0; beta_b(1) = 0
      phip_b(1) = 0; phi_b(1) = 0; bvco_b(1) = 0
      buco_b(1) = 0

      DO nsval = 2, ns_b
         READ(iunit, iostat=ierr, err=100) iota_b(nsval),               &
         pres_b(nsval), beta_b(nsval), phip_b(nsval), phi_b(nsval),     &
         bvco_b(nsval), buco_b(nsval)
      END DO

      READ(iunit, iostat=ierr, err=100) mboz_b, nboz_b, mnboz_b, jsize
      READ(iunit, iostat=js) version, lasym_b

      ALLOCATE (bmnc_b(mnboz_b,ns_b), rmnc_b(mnboz_b,ns_b),             &
        zmns_b(mnboz_b,ns_b), pmns_b(mnboz_b,ns_b),                     &
        gmnc_b(mnboz_b,ns_b), ixm_b(mnboz_b), ixn_b(mnboz_b), stat=ierr)
      IF (ierr .ne. 0) THEN
         PRINT *,' Allocation error in read_boozer_file'
         RETURN
      END IF

      idx_b = 0
      ixm_b = 0
      rmnc_b = 0; zmns_b = 0; pmns_b = 0; bmnc_b = 0; gmnc_b = 0

      IF (lasym_b) THEN
      ALLOCATE (bmns_b(mnboz_b,ns_b), rmns_b(mnboz_b,ns_b),             &
        zmnc_b(mnboz_b,ns_b), pmnc_b(mnboz_b,ns_b),                     &
        gmns_b(mnboz_b,ns_b), stat=ierr)
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

        READ(iunit,iostat=ierr,err=100, END=200) bmnc_b(:mnboz_b,nsval),   &
             rmnc_b(:mnboz_b,nsval), zmns_b(:mnboz_b,nsval),               &
             pmns_b(:mnboz_b,nsval), gmnc_b(:mnboz_b,nsval)

        IF (.not.lasym_b) CYCLE

        READ(iunit,iostat=ierr,err=100, END=200) bmns_b(:mnboz_b,nsval),   &
             rmns_b(:mnboz_b,nsval), zmnc_b(:mnboz_b,nsval),            &
             pmnc_b(:mnboz_b,nsval), gmns_b(:mnboz_b,nsval)
      END DO

 100  CONTINUE
      IF (ierr .gt. 0) THEN
         PRINT *,' Error reading in subroutine read_boozer_file:',      &
                 ' ierr = ', ierr
      END IF
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

      IF (ALLOCATED(iota_b)) DEALLOCATE (iota_b, pres_b, beta_b,        &
          phip_b, phi_b, bvco_b, buco_b, idx_b, stat = istat)

      IF (ALLOCATED(bmnc_b)) DEALLOCATE (bmnc_b, rmnc_b,                &
         zmns_b, pmns_b, gmnc_b, ixm_b, ixn_b, stat = istat)

      IF (ALLOCATED(bmns_b)) DEALLOCATE (bmns_b, rmns_b,                &
         zmnc_b, pmnc_b, gmns_b, stat = istat)

      END SUBROUTINE read_boozer_deallocate


#if defined(NETCDF)
      SUBROUTINE unpack_cdf (nbooz, var_name, array2d)
      USE stel_kinds, ONLY: rprec
      USE ezcdf
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: nbooz
      INTEGER :: nsval, icount
      REAL(rprec), DIMENSION(:,:), INTENT(out) :: array2d
      CHARACTER(LEN=*), INTENT(in) :: var_name
!-----------------------------------------------
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
#endif

      END MODULE read_boozer_mod
