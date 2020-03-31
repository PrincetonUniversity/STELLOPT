      SUBROUTINE write_boozmn (extension)
      USE booz_params
      USE booz_persistent, ONLY: xmb, xnb 
      USE read_wout_mod, ONLY: rmax_surf, rmin_surf, betaxis, aspect,
     1                         lasym
!DEC$ IF DEFINED (NETCDF)
      USE ezcdf
C-----------------------------------------------
C   L O C A L   P A R A M E T E R S
C-----------------------------------------------
      USE read_boozer_mod, ONLY: vn_nfp, vn_ns, vn_aspect,
     1 vn_rmax, vn_rmin, vn_betaxis, vn_mboz, vn_nboz, vn_mnboz,
     2 vn_version, vn_iota, vn_pres, vn_beta, vn_phip, vn_phi,
     3 vn_bvco, vn_buco, vn_ixm, vn_ixn, vn_bmnc, vn_rmnc,
     4 vn_zmns, vn_pmns, vn_gmnc, vn_lasym, vn_bmns, vn_rmns,
     5 vn_zmnc, vn_pmnc, vn_gmns, vn_jlist
!DEC$ ELSE
      USE safe_open_mod
!DEC$ ENDIF
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER(LEN=*) :: extension
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      CHARACTER(LEN=*), PARAMETER :: version =
     1   "Boozer Transformation Code Version 2.0"
!-----------------------------------------------
!DEC$ IF DEFINED (NETCDF)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: nbooz, ierr
      CHARACTER(LEN=*), PARAMETER, DIMENSION(1) ::
     1             r1dim = (/'radius'/), mn1dim = (/'mn_mode'/),
     2             j1dim = (/'comput_surfs'/)
      CHARACTER(LEN=*), DIMENSION(2), PARAMETER :: 
     1             r2dim = (/'mn_modes','pack_rad'/)
C-----------------------------------------------
! Open cdf File
      CALL cdf_open(nbooz,'boozmn_' // TRIM(extension) // '.nc',
     1             'w', ierr)
      IF (ierr .ne. 0) THEN
         PRINT *,' Error opening boozmn .nc file'
         RETURN
      END IF

! Define Variables
! Scalars
      CALL cdf_define(nbooz, vn_nfp, nfp)
      CALL cdf_define(nbooz, vn_ns, ns)
      CALL cdf_define(nbooz, vn_aspect, aspect)
      CALL cdf_define(nbooz, vn_rmax, rmax_surf)
      CALL cdf_define(nbooz, vn_rmin, rmin_surf)
      CALL cdf_define(nbooz, vn_betaxis, betaxis)
      CALL cdf_define(nbooz, vn_mboz, mboz)
      CALL cdf_define(nbooz, vn_nboz, nboz)
      CALL cdf_define(nbooz, vn_mnboz, mnboz)
      CALL cdf_define(nbooz, vn_version, version)
      CALL cdf_define(nbooz, vn_lasym, lasym_b)
! 1D Arrays
      CALL cdf_define(nbooz, vn_iota, hiota, dimname=r1dim)
      CALL cdf_define(nbooz, vn_pres, pres, dimname=r1dim)
      CALL cdf_define(nbooz, vn_beta, beta_vol, dimname=r1dim)
      CALL cdf_define(nbooz, vn_phip, phip, dimname=r1dim)
      CALL cdf_define(nbooz, vn_phi,  phi, dimname=r1dim)
      CALL cdf_define(nbooz, vn_bvco, bvco, dimname=r1dim)
      CALL cdf_define(nbooz, vn_buco, buco, dimname=r1dim)
      CALL cdf_define(nbooz, vn_jlist, jlist, dimname=j1dim)
      CALL cdf_define(nbooz, vn_ixm, NINT(xmb), dimname=mn1dim)
      CALL cdf_define(nbooz, vn_ixn, NINT(xnb), dimname=mn1dim)
! 2D Arrays
      CALL cdf_define(nbooz, vn_bmnc, bmncb, dimname=r2dim)
      CALL cdf_define(nbooz, vn_rmnc, rmncb, dimname=r2dim)
      CALL cdf_define(nbooz, vn_zmns, zmnsb, dimname=r2dim)
      CALL cdf_define(nbooz, vn_pmns, pmnsb, dimname=r2dim)
      CALL cdf_define(nbooz, vn_gmnc, gmncb, dimname=r2dim)
      IF (lasym_b) THEN
      CALL cdf_define(nbooz, vn_bmns, bmnsb, dimname=r2dim)
      CALL cdf_define(nbooz, vn_rmns, rmnsb, dimname=r2dim)
      CALL cdf_define(nbooz, vn_zmnc, zmncb, dimname=r2dim)
      CALL cdf_define(nbooz, vn_pmnc, pmncb, dimname=r2dim)
      CALL cdf_define(nbooz, vn_gmns, gmnsb, dimname=r2dim)
      ENDIF

! Write out scalars
      CALL cdf_write(nbooz, vn_nfp, nfp)
      CALL cdf_write(nbooz, vn_ns, ns)
      CALL cdf_write(nbooz, vn_aspect, aspect)
      CALL cdf_write(nbooz, vn_rmax, rmax_surf)
      CALL cdf_write(nbooz, vn_rmin, rmin_surf)
      CALL cdf_write(nbooz, vn_betaxis, betaxis)
      CALL cdf_write(nbooz, vn_mboz, mboz)
      CALL cdf_write(nbooz, vn_nboz, nboz)
      CALL cdf_write(nbooz, vn_mnboz, mnboz)
      CALL cdf_write(nbooz, vn_version, version)

!  1D arrays 
      hiota(1) = 0; pres(1) = 0; beta_vol(1) = 0
      phip(1) = 0; phi(1) = 0; bvco(1) = 0; buco(1) = 0

      CALL cdf_write(nbooz, vn_iota, hiota)
      CALL cdf_write(nbooz, vn_pres, pres)
      CALL cdf_write(nbooz, vn_beta, beta_vol)
      CALL cdf_write(nbooz, vn_phip, phip)
      CALL cdf_write(nbooz, vn_phi,  phi)
      CALL cdf_write(nbooz, vn_bvco, bvco)
      CALL cdf_write(nbooz, vn_buco, buco)
      CALL cdf_write(nbooz, vn_jlist, jlist)
      CALL cdf_write(nbooz, vn_ixm, NINT(xmb))
      CALL cdf_write(nbooz, vn_ixn, NINT(xnb))

!  Write packed 2D arrays
      CALL cdf_write(nbooz, vn_bmnc, bmncb)
      CALL cdf_write(nbooz, vn_rmnc, rmncb)
      CALL cdf_write(nbooz, vn_zmns, zmnsb)
      CALL cdf_write(nbooz, vn_pmns, pmnsb)
      CALL cdf_write(nbooz, vn_gmnc, gmncb)
      CALL cdf_write(nbooz, vn_lasym, lasym_b)
      IF (lasym_b) THEN
      CALL cdf_write(nbooz, vn_bmns, bmnsb)
      CALL cdf_write(nbooz, vn_rmns, rmnsb)
      CALL cdf_write(nbooz, vn_zmnc, zmncb)
      CALL cdf_write(nbooz, vn_pmnc, pmncb)
      CALL cdf_write(nbooz, vn_gmns, gmnsb)
      ENDIF


! Close cdf File
      CALL cdf_close(nbooz, ierr)

!DEC$ ELSE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, iunit, istat, js
C-----------------------------------------------
      iunit = unit_booz
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
     1   nfp, ns, aspect, rmax_surf, rmin_surf, betaxis

      DO js = 2, ns
         WRITE(iunit, iostat=istat, err=100) hiota(js), pres(js),
     1   beta_vol(js), phip(js), phi(js), bvco(js), buco(js)
      END DO

!     SPH 070909: ADDED lasym to dump
      WRITE (iunit, iostat=istat, err=100) mboz, nboz, mnboz, jsize
      WRITE (iunit, iostat=istat, err=100) version, lasym

      WRITE (iunit, iostat=istat, err=100) NINT(xnb(:mnboz)), 
     1                                     NINT(xmb(:mnboz))

!
!     Write packed (in radial coordinate) 2D arrays
!
      DO i = 1, jsize
        js = jlist(i)
        IF (js.le.0 .or. js.gt.ns) CYCLE
        WRITE (iunit, iostat=istat, err=100) js
        WRITE (iunit, iostat=istat, err=100) bmncb(:mnboz,i),
     1       rmncb(:mnboz,i), zmnsb(:mnboz,i), pmnsb(:mnboz,i),
     2       gmncb(:mnboz,i)
!SPH070909: WRITE OUT ASYMMETRIC PARTS
        IF (lasym) THEN
        WRITE (iunit, iostat=istat, err=100) bmncb(:mnboz,i),
     1       rmnsb(:mnboz,i), zmncb(:mnboz,i), pmncb(:mnboz,i),
     2       gmnsb(:mnboz,i)
        ENDIF	
      END DO


 100  CONTINUE
      IF (istat .gt. 0)
     1    PRINT *,' Error writing in subroutine write_boozmn:',
     2            ' istat = ', istat

      CLOSE(iunit)

!DEC$ ENDIF
      END SUBROUTINE write_boozmn
