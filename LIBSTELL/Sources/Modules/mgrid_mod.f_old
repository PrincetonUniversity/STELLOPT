      MODULE mgrid_mod
      USE stel_kinds
      USE vmec_input, ONLY: nbfld, nflxs, lfreeb, lrecon
      USE vsvd0, ONLY: nigroup, nparts, npfcoil, nbcoilsp, nfloops,
     1                 nbctotp
      IMPLICIT NONE
       
      LOGICAL :: lnverror=.true.
      INTEGER, PARAMETER :: nlimset = 2       !number of different limiters
      CHARACTER(LEN=*), PARAMETER :: 
     1   vn_br0 = 'br', vn_bp0 = 'bp', vn_bz0 = 'bz',
     3   vn_ir = 'ir', vn_jz = 'jz',
     4   vn_kp = 'kp', vn_nfp = 'nfp',
     5   vn_rmin='rmin', vn_rmax='rmax', vn_zmin='zmin', 
     6   vn_zmax='zmax', vn_coilgrp='coil_group' 
      CHARACTER(LEN=*), PARAMETER ::
     1  vn_nextcur = 'nextcur',  vn_mgmode='mgrid_mode', 
     2  vn_coilcur = 'raw_coil_cur',
     3  vn_flp = 'nobser', vn_nobd = 'nobd', vn_nbset = 'nbsets', 
     4  vn_nbfld = 'nbfld',
     2  ln_flp = 'flux loops', ln_nobd = 'Connected flux loops',
     3  ln_nbset = 'B-coil loops', ln_next = 'External currents',
     4  ln_nbfld = 'B-coil measurements'
 
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
!     nr0b, np0b, nz0b:
!              :  grid dimensions of magnetic field in mgrid file
!     nbvac    :  total number of grid points (nr0b*np0b*nz0b) in mgrid file
!     bvac(:,1):  br (radial component of external magnetic field)
!     bvac(:,2):  bp (toroidal component)
!     bvac(:,3) = bz (z-component)
!     rminb, rmaxb : min (max) radial dimension of grid in mgrid
!     zminb, zmaxb : min (max) vertical dimension of grid in mgrid
!
!     nextcur:         no. of EXTERNAL current groups (eg., TF, PF, helical)
!     raw_coil_current  array of raw currents for each coil group
!     mgrid_mode     = 'S', scaled mode; = 'R', raw mode
!     curlabel:   array of labels describing each current group
!                     included in green''s FUNCTION BFIELD response
!
!               - - - - - - - - - - - - - - - - - -
!               FOR DIAGNOSTICS AND DATA ANALYSIS
!               (HERE,COILS ARE FOR MEASURING FIELDS, FLUXES)
!    iconnect:   two-dimensional array describing electrical
!     needflx:    =NEEDIT, loop required for flux match
!                    >=ISYMCOIL, loop required but flux computed by
!                               invoking symmetry in Z
!                    =IDONTNEED, loop not required for flux match
!    needbfld:    =NEEDIT, loop required for B-field match
!                    =ISAMECOIL, loop at same position as previous loop
!                    >=ISYMCOIL, loop required but B-field computed
!                               by invoking symmetry in Z
!                    =IDONTNEED, loop not required in B-field match
!      dsiext:    connected flux loop signals due to EXTERNAL coils
!      plflux:    array of measured (inferred) plasma contrib. to flux loops
!      plbfld:    array of measured (inferred) plasma contrib. to B-loops
!                      connection of up to four flux loops. Specifies
!                     the sign and flux loop number of (up to) four
!                     connected individual loops (indexing based on
!                     xobser,zobser arrays).
!        nobd:   number of connected flux loop measurements
!      nobser:   number of individual flux loop positions
!      nbsets:   number of B-coil sets defined in mgrid file
!  nbcoils(n):   number of bfield coils in each set defined in mgrid file
!    nbcoilsn:   total number of bfield coils defined in mgrid file
!      nbfldn:   total number of EXTERNAL bfield measurements used in matching
!  bloopnames:   array of labels describing b-field sets
!    dsilabel:   array of labels describing connected flux loops
!      xobser:   array of flux loop R-positions
!      zobser:   array of flux loop Z-positions
! rbcoil(m,n):   R position of the m-th coil in the n-th set from mgrid file
! zbcoil(m,n):   Z position of the m-th coil in the n-th set from mgrid file
! abcoil(m,n):   orientation (surface normal wrt R axis; in radians)
!
      INTEGER :: nr0b, np0b, nfper0, nz0b
      INTEGER :: nobd, nobser, nextcur, nbfldn, nbsets, nbcoilsn
      INTEGER :: nbvac, nbcoil_max, nlim, nlim_max, nsets,
     1           nrgrid, nzgrid
      INTEGER, DIMENSION(:), ALLOCATABLE :: needflx, nbcoils
      INTEGER, DIMENSION(:), ALLOCATABLE :: limitr, nsetsn
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: iconnect, needbfld
      REAL(rprec) :: rminb, zminb, rmaxb, zmaxb, delrb, delzb
      REAL(rprec) ::rx1, rx2, zy1, zy2, condif
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE, TARGET :: bvac
      REAL(rprec), DIMENSION(:,:,:), POINTER :: brvac, bzvac, bpvac
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: unpsiext, 
     1   plbfld, rbcoil, zbcoil, abcoil, bcoil, rbcoilsqr
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: raw_coil_current
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: xobser, zobser,
     1   xobsqr, dsiext, psiext, plflux, b_chi
      CHARACTER(LEN=300) :: mgrid_path
      CHARACTER(LEN=300) :: mgrid_path_old = " "
      CHARACTER(LEN=30), DIMENSION(:), ALLOCATABLE :: curlabel
      CHARACTER(LEN=15), DIMENSION(:), ALLOCATABLE :: 
     1                                           dsilabel, bloopnames
      CHARACTER(LEN=30) :: tokid
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: dbcoil, pfcspec
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: 
     1    rlim, zlim, reslim, seplim
      CHARACTER(LEN=1) :: mgrid_mode

!DEC$ IF DEFINED (NETCDF)
      PRIVATE :: read_mgrid_bin, read_mgrid_nc
!DEC$ ELSE
      PRIVATE :: read_mgrid_bin
!DEC$ ENDIF
      
      CONTAINS

      SUBROUTINE read_mgrid (mgrid_file, extcur, nv, nfp, lscreen, 
     1                       ier_flag)
      USE system_mod
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   V a r i a b l e s
C-----------------------------------------------
!     
!     mgrid_file:     full path to mgrid file
!     lscreen   :     logical controlling output to screen
!     ier_flag  ;     error flag returned to caller
!     extcur(n)    :  external current multiplier for bfield(n) components
!
      INTEGER, INTENT(out) :: ier_flag
      INTEGER, INTENT(in)  :: nv, nfp
      LOGICAL, INTENT(in)  :: lscreen
      REAL(rprec), INTENT(in) :: extcur(:)
      CHARACTER(len=*), INTENT(in) :: mgrid_file
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
!DEC$ IF DEFINED (VMS)
      CHARACTER(LEN=*), PARAMETER :: mgrid_defarea='vmec$:[makegrid]'
!DEC$ ELSE
      CHARACTER(LEN=*), PARAMETER :: mgrid_defarea='$HOME/vmec/MAKEGRID'
!DEC$ ENDIF
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
!
!     lgrid_exist  :   logical set if mgrid file is found in given path
!
      INTEGER :: istat, ii
      CHARACTER(LEN=200) :: home_dir
      LOGICAL :: lgrid_exist, lfind
C-----------------------------------------------
      mgrid_path(:) = ' '
      mgrid_path = TRIM(mgrid_file)

      IF ((mgrid_path .eq. TRIM(mgrid_path_old)) .and. 
     1     ALLOCATED(curlabel)) THEN
         PRINT *,' mgrid file previously parsed!'
         RETURN
      END IF
      
      INQUIRE (file=mgrid_path,exist=lgrid_exist,iostat=istat)
      IF (istat.ne.0 .or. .not.lgrid_exist) THEN
          IF (lscreen) PRINT *,' MGRID FILE NOT FOUND IN SPECIFIED ',
     1       'PATH: SEARCHING DEFAULT AREA'
          ii = INDEX(mgrid_file,'/',back=.true.)
          istat = INDEX(mgrid_defarea, '$HOME')
          IF (istat .ne. 0) THEN
             CALL getenv('HOME', home_dir)
             IF (istat .gt. 1) THEN
                home_dir = mgrid_defarea(1:istat-1) // TRIM(home_dir)
     1                   // mgrid_defarea(istat+5:)
             ELSE
                home_dir = TRIM(home_dir) // mgrid_defarea(istat+5:)
             END IF
          ELSE
             home_dir = mgrid_defarea
          END IF
          mgrid_path = TRIM(home_dir) // mgrid_file(ii+1:)
          INQUIRE (file=mgrid_path,exist=lgrid_exist,iostat=istat)
      END IF

      mgrid_path_old = mgrid_path

      ier_flag = 0

      IF (lgrid_exist) THEN
         IF (lscreen) PRINT '(2x,2a)',
     1     'Opening vacuum field file: ', TRIM(mgrid_file)
!
!        Parse mgrid file name, look for .nc extension (netcdf format)
! 
         ii = LEN_TRIM(TRIM(mgrid_path)) - 2
         lfind = (mgrid_path(ii:ii+2) == '.nc')
         IF (lfind) THEN
!DEC$ IF DEFINED (NETCDF)
            CALL read_mgrid_nc (mgrid_path, extcur, nv, nfp, 
     1                          ier_flag, lscreen)
!DEC$ ELSE
            lgrid_exist = .false.
!DEC$ ENDIF
         ELSE
            CALL read_mgrid_bin (mgrid_path, extcur, nv, nfp,
     1                          ier_flag, lscreen)
         END IF

         IF (np0b .ne. nv) THEN
            IF (lnverror) PRINT *,' NZETA=',nv,
     1      ' NOT EQUAL TO NP0B=',np0b,' IN MGRID FILE'
            ier_flag = 9
         ELSE IF (nfper0.ne.nfp) THEN
            PRINT *,' NFP(READ in) = ',nfp,' DOES NOT AGREE WITH ',
     1      'NFPER (in vacuum field file) = ',nfper0
            ier_flag = 9
         END IF

      END IF
      
      IF (ier_flag .ne. 0) RETURN

      IF (.not.lgrid_exist .or. ier_flag.ne.0) THEN
         lfreeb = .false.
         lrecon = .false.
         IF (lscreen) THEN
            PRINT *, ' Error opening/reading mgrid file in dir: ',
     1                TRIM(home_dir)
            PRINT *, ' User must supply vacuum bfield in mgrid to ',
     1                'run vmec in free-boundary mode!'
            PRINT *, ' Proceeding to run vmec in',
     1               ' fixed boundary mode'
         END IF
      END IF

      END SUBROUTINE read_mgrid

      
      SUBROUTINE read_mgrid_bin (filename, extcur, nv, nfp, ier_flag, 
     1                           lscreen)
      USE safe_open_mod
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y  A r g u m e n t s
C-----------------------------------------------
!
!     lstyle_2000  :   logical controlling ordering of magnetic field components read in
!
      INTEGER, INTENT(in)  :: nv, nfp
      CHARACTER(LEN=*), INTENT(in) :: filename
      REAL(rprec), INTENT(in) :: extcur(:)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: 
     1           brtemp, bztemp, bptemp
      INTEGER :: ier_flag, iunit = 50
      INTEGER :: istat, ig, i, j, n, n1, m, nsets_max, k 
      LOGICAL :: lscreen, lstyle_2000
C-----------------------------------------------

      CALL safe_open(iunit, istat, filename, 'old', 'unformatted')
      IF (istat .ne. 0) THEN
         ier_flag = 9
         RETURN
      END IF

      READ (iunit,iostat=istat) nr0b, nz0b, np0b, nfper0, nextcur
      IF (istat .ne. 0) ier_flag = 9

      IF (nfper0.ne.nfp .or. np0b.ne.nv) RETURN

      lstyle_2000 = (nextcur < 0)
      nextcur = ABS(nextcur)
      READ(iunit,iostat=istat) rminb, zminb, rmaxb, zmaxb
      IF (istat .ne. 0) ier_flag = 9

      IF (nextcur .eq. 0) THEN
        PRINT *,' NEXTCUR = 0 IN READING MGRID FILE'
        ier_flag = 9
      ELSE IF (nextcur .gt. nigroup) THEN
        PRINT *,' NEXTCUR > NIGROUP IN MGRID FILE'
        ier_flag = 9
      END IF

      IF (ier_flag .ne. 0) RETURN

      ALLOCATE (curlabel(5*(nextcur/5+1)), stat=istat)    !MIN of 5 for printing
      curlabel = " "
      READ(iunit,iostat=istat) (curlabel(n),n=1,nextcur)
      IF (istat .ne. 0) THEN
         PRINT *,' reading mgrid file failed (curlabel)'
         ier_flag = 9
         RETURN
      END IF

!
!     NOTE: ADD UP CONTRIBUTION TO BVAC DIRECTLY FOR ALL EXTERNAL CURRENT GROUPS

      nbvac = nr0b*nz0b*np0b
      IF (.NOT. ALLOCATED(bvac)) THEN
         ALLOCATE (bvac(nbvac,3))
      ELSE IF (SIZE(bvac,1) .ne. nbvac) THEN
         DEALLOCATE (bvac);  ALLOCATE(bvac(nbvac,3))
      END IF

      ALLOCATE (brtemp(nr0b,nz0b,np0b), bptemp(nr0b,nz0b,np0b),
     1          bztemp(nr0b,nz0b,np0b), stat=istat)

      IF (istat .ne. 0) THEN
        PRINT *,' allocation for b-vector storage failed'
        ier_flag = 9
        RETURN
      END IF

      bvac = 0

      DO ig = 1,nextcur
         IF (lstyle_2000) THEN
            READ(iunit, iostat=istat) brtemp, bptemp, bztemp
         ELSE
            READ(iunit, iostat=istat) (((brtemp(i,j,k), bztemp(i,j,k),
     1                                   bptemp(i,j,k), i= 1,nr0b), 
     2                                   j=1,nz0b), k=1,np0b)
         END IF
!
!        STORE SUMMED BFIELD (OVER COIL GROUPS) IN BVAC
!         
         CALL sum_bfield(bvac(1,1), brtemp, extcur(ig), nbvac)
         CALL sum_bfield(bvac(1,2), bptemp, extcur(ig), nbvac)
         CALL sum_bfield(bvac(1,3), bztemp, extcur(ig), nbvac)

      END DO

      DEALLOCATE (brtemp, bztemp, bptemp)

      IF (istat .ne. 0) THEN
         ier_flag = 9
         RETURN
      END IF

      IF (lstyle_2000) THEN
         READ (iunit, iostat=istat) mgrid_mode
         IF (istat .eq. 0) THEN
            ALLOCATE (raw_coil_current(nextcur))
            READ (iunit, iostat=istat) raw_coil_current(1:nextcur)
            IF (istat .ne. 0) mgrid_mode = 'N'
         END IF
      ELSE
         mgrid_mode = 'N'         !Old-style, no mode info
      END IF

!
!     READ IN EXTERNAL POLOIDAL FLUX, FIELD MEASURMENT
!     LOOP COORDINATES AND LABELS
!
      READ(iunit,iostat=istat) nobser, nobd, nbsets
      IF (istat.ne.0) THEN
         nobser = 0
         nobd   = 0
         nbsets = 0
         IF (lscreen) PRINT *,' No observation data in mgrid data'
         GOTO 900
      END IF

      nbfldn = SUM(nbfld(:nbsets))
      ALLOCATE (nbcoils(nbsets), stat=istat)
      READ(iunit) (nbcoils(n),n=1,nbsets)

      nbcoil_max = MAXVAL(nbcoils(:nbsets))

      ALLOCATE (xobser(nobser), zobser(nobser), dsilabel(nobd),
     1       iconnect(4,nobser+nobd), unpsiext(nobser,nextcur),
     2       xobsqr(nobser), needflx(nobser), plflux(nobser+nobd),
     3       dsiext(nobd), psiext(nobser), bloopnames(nbsets),
     4       needbfld(nbcoil_max,nbsets), plbfld(nbcoil_max,nbsets),
     5       rbcoil(nbcoil_max,nbsets), zbcoil(nbcoil_max,nbsets),
     6       abcoil(nbcoil_max,nbsets), bcoil(nbcoil_max,nbsets),
     7       rbcoilsqr(nbcoil_max,nbsets), b_chi(nbsets),
     8       dbcoil(nbcoil_max,nbsets,nextcur), stat = istat)
      IF (istat .ne. 0) THEN
          IF (lscreen)
     1       PRINT *,' allocation error for xobser: istat = ',istat
          ier_flag = 9
          RETURN
      END IF

      IF (nobser .gt. nfloops) THEN
         PRINT *, 'NOBSER>NFLOOPS'
         ier_flag = 9
      END IF
      IF (nobd .gt. nfloops) THEN
         PRINT *, 'NOBD>NFLOOPS'
         ier_flag = 9
      END IF
      IF (nflxs .gt. nfloops) THEN
         PRINT *, 'NFLXS>NFLOOPS'
         ier_flag = 9
      END IF
      IF (nbfldn .gt. nbctotp) THEN
         PRINT *, 'NBFLDN>NBCTOTP'
         ier_flag = 9
      END IF
      IF (nbcoil_max .gt. nbcoilsp) THEN
         PRINT *, 'NBCOIL_max>NBCOILSP'
         ier_flag = 9
      END IF

      IF (ier_flag .ne. 0) RETURN

      IF (nobser+nobd .gt. 0) iconnect(:4,:nobser+nobd) = 0

      READ(iunit) (xobser(n), zobser(n),n=1,nobser)
      READ(iunit) (dsilabel(n),n=1,nobd)
      READ(iunit) ((iconnect(j,n),j=1,4),n=1,nobd)

      IF (nbcoil_max.gt.0 .and. nbsets.gt.0) THEN
         rbcoil(:nbcoil_max,:nbsets) = 0
         zbcoil(:nbcoil_max,:nbsets) = 0
         abcoil(:nbcoil_max,:nbsets) = 0

         DO n=1,nbsets
           IF (nbcoils(n).gt.0) THEN
           READ(iunit) n1,bloopnames(n1)
           READ(iunit)(rbcoil(m,n),zbcoil(m,n),abcoil(m,n),
     1             m=1,nbcoils(n))
           ENDIF
         ENDDO

         dbcoil(:nbcoil_max,:nbsets,:nextcur) = 0
      END IF
      DO ig = 1,nextcur
        !un-connected coil fluxes
         READ(iunit) (unpsiext(n,ig),n=1,nobser)
         DO n = 1,nbsets
            READ(iunit) (dbcoil(m,n,ig),m=1,nbcoils(n))
         ENDDO
      ENDDO

!
!     READ LIMITER & PROUT PLOTTING SPECS
!
      ALLOCATE (limitr(nlimset), nsetsn(nigroup))

      READ (iunit,iostat=istat) nlim,(limitr(i),i=1,nlim)
      IF (istat .ne. 0)then
        nlim = 0
        IF (lscreen) PRINT *,' No limiter data in mgrid file'
        GOTO 900
      END IF

      nlim_max = MAXVAL(limitr)

      IF (nlim .gt. nlimset) THEN
         PRINT *, 'nlim>nlimset'
         ier_flag = 9
         RETURN
      END IF

      ALLOCATE( rlim(nlim_max,nlim),   zlim(nlim_max,nlim),
     1          reslim(nlim_max,nlim) ,seplim(nlim_max,nlim),
     2          stat=istat)
      IF (istat .ne. 0) THEN
         PRINT *, 'rlim istat!=0'
         ier_flag = 9
         RETURN
      END IF

      READ(iunit, iostat=istat)
     1   ((rlim(i,j),zlim(i,j),i=1,limitr(j)),j=1,nlim)
      READ(iunit, iostat=istat) nsets,(nsetsn(i), i=1,nsets)

      IF (nsets .gt. nigroup) THEN
         PRINT *, 'nsets>nigroup'
         ier_flag = 9
         RETURN
      ELSE IF (istat .ne. 0) THEN
         ier_flag = 9
         RETURN
      END IF

      nsets_max = MAXVAL(nsetsn)

      IF (nsets_max .gt. npfcoil) THEN
         PRINT *, 'nsetsn>npfcoil'
         ier_flag = 9
         RETURN
      END IF

      ALLOCATE (pfcspec(nparts,nsets_max,nsets), stat=istat)

!     NOTE TO RMW: SHOULD READ IN NPARTS HERE (PUT INTO MGRID FILE)

      READ(iunit, iostat=istat) (((pfcspec(i,j,k),i=1,nparts),
     1        j=1,nsetsn(k)), k=1,nsets)

      DEALLOCATE (limitr, nsetsn)

      READ(iunit, iostat=istat) rx1,rx2,zy1,zy2,condif,
     1  nrgrid,nzgrid,tokid

      IF (istat .ne. 0) THEN
         ier_flag = 9
         RETURN
      END IF

      IF (nobser .gt. 0) xobsqr(:nobser) = SQRT(xobser(:nobser))
!
!       PARTITION MGRID B-LOOPS INTO SETS
!
      nbcoilsn = SUM(nbcoils(:nbsets))

      DO n = 1,nbsets
        rbcoilsqr(:nbcoils(n),n) = SQRT(rbcoil(:nbcoils(n),n))
      ENDDO

 900  CONTINUE

      CLOSE (iunit)

      delrb = (rmaxb-rminb)/(nr0b-1)
      delzb = (zmaxb-zminb)/(nz0b-1)

!
!     SUM UP CONTRIBUTIONS FROM INDIVIDUAL COIL GROUPS
!
      IF (lfreeb) THEN
         IF (nobser .gt. 0) psiext(:nobser) = 0
         IF (nbcoil_max.gt.0 .and. nbsets.gt.0)
     1       bcoil(:nbcoil_max, :nbsets) = 0

         DO ig = 1,nextcur
            IF (nobser .gt. 0)
     1      psiext(:nobser) = psiext(:nobser) +
     2                        extcur(ig)*unpsiext(:nobser,ig)
            DO n=1,nbsets
               n1 = nbcoils(n)
               bcoil(:n1,n) = bcoil(:n1,n) +
     1                        extcur(ig)*dbcoil(:n1,n,ig)
            ENDDO
         ENDDO
      ENDIF                   !!IF LFREEB

      END SUBROUTINE read_mgrid_bin

!DEC$ IF DEFINED (NETCDF)
      SUBROUTINE read_mgrid_nc (filename, extcur, nv, nfp,
     1                          ier_flag, lscreen)
      USE ezcdf
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y  A r g u m e n t s
C-----------------------------------------------
      CHARACTER(LEN=*), INTENT(in) :: filename
      INTEGER, INTENT(in)     :: nv, nfp
      REAL(rprec), INTENT(in) :: extcur(:)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: 
     1                       brtemp, bztemp, bptemp
      INTEGER :: ier_flag, ngrid
      INTEGER :: istat, ig
      LOGICAL :: lscreen
      INTEGER, DIMENSION(3)   :: dimlens
      CHARACTER(LEN=100) :: temp
      CHARACTER(LEN=100), ALLOCATABLE, DIMENSION(:) :: 
     1                 vn_br, vn_bp, vn_bz
C-----------------------------------------------
      call cdf_open(ngrid, filename,'r', istat)
      IF (istat .ne. 0) THEN
         ier_flag = 9
         RETURN
      END IF

!
!     READ IN DATA
!
      CALL cdf_read(ngrid, vn_ir, nr0b)
      CALL cdf_read(ngrid, vn_jz, nz0b)
      CALL cdf_read(ngrid, vn_kp, np0b)
      CALL cdf_read(ngrid, vn_nfp, nfper0)

      IF (nfper0.ne.nfp .or. np0b.ne.nv) RETURN

      CALL cdf_read(ngrid, vn_nextcur, nextcur)

      IF (nextcur .eq. 0) THEN
        PRINT *,' NEXTCUR = 0 IN READING MGRID FILE'
        ier_flag = 9
        RETURN
      ELSE IF (nextcur .gt. nigroup) THEN
        PRINT *,' NEXTCUR > NIGROUP IN MGRID FILE'
        ier_flag = 9
        RETURN
      END IF

      CALL cdf_read(ngrid, vn_rmin, rminb)
      CALL cdf_read(ngrid, vn_zmin, zminb)
      CALL cdf_read(ngrid, vn_rmax, rmaxb)
      CALL cdf_read(ngrid, vn_zmax, zmaxb)

      delrb = (rmaxb-rminb)/(nr0b-1)
      delzb = (zmaxb-zminb)/(nz0b-1)

      CALL cdf_inquire(ngrid, vn_coilgrp, dimlens)
	IF (.NOT. ALLOCATED(curlabel)) THEN
          ALLOCATE (curlabel(nextcur), stat=istat)
	ELSE IF (SIZE(curlabel) .ne. nextcur) THEN
	    DEALLOCATE (curlabel)
          ALLOCATE (curlabel(nextcur), stat=istat)
	END IF
!THIS IS A GLITCH WITH cdf_read: must distinguish 1D char array from multi-D
      IF (nextcur .eq. 1) THEN
         IF (istat .eq. 0) 
     1     CALL cdf_read(ngrid, vn_coilgrp, curlabel(1))
      ELSE IF (istat .eq. 0) THEN
           CALL cdf_read(ngrid, vn_coilgrp, curlabel(1:nextcur))
      END IF

	IF (istat .ne. 0) STOP 'Error allocating CURLABEL in mgrid_mod'
!
!     READ 3D Br, Bp, Bz ARRAYS FOR EACH COIL GROUP
!
      ALLOCATE (vn_br(nextcur), vn_bz(nextcur), vn_bp(nextcur), 
     1          stat=istat)
	IF (istat .ne. 0) STOP 'Error allocating vn_bX in mgrid_mod'

      nbvac = nr0b*nz0b*np0b
      IF (.NOT. ALLOCATED(bvac)) THEN
         ALLOCATE (bvac(nbvac,3), stat=istat)
      ELSE IF (SIZE(bvac,1) .ne. nbvac) THEN
         DEALLOCATE (bvac);  ALLOCATE(bvac(nbvac,3), stat=istat)
      END IF
      IF (istat .ne. 0) STOP 'Error allocating bvac in mgrid_mod'

      bvac = 0

      DO ig = 1, nextcur
         WRITE (temp, '(a,i3.3)') "_",ig
         vn_br(ig) = vn_br0 // temp
         vn_bp(ig) = vn_bp0 // temp
         vn_bz(ig) = vn_bz0 // temp
         CALL cdf_inquire(ngrid, vn_br(ig), dimlens)
         IF (.NOT. ALLOCATED(brtemp)) THEN
            ALLOCATE (brtemp(dimlens(1),dimlens(2),dimlens(3)),
     1                bptemp(dimlens(1),dimlens(2),dimlens(3)),
     2                bztemp(dimlens(1),dimlens(2),dimlens(3)),
     3                stat=istat)
            IF (istat .ne. 0)STOP 'Error allocating bXtemp in mgrid_mod'
         END IF
         CALL cdf_read(ngrid, vn_br(ig), brtemp)
         CALL cdf_read(ngrid, vn_bp(ig), bptemp)
         CALL cdf_read(ngrid, vn_bz(ig), bztemp)

!
!        STORE SUMMED BFIELD (OVER COIL GROUPS) IN BVAC
!         
         CALL sum_bfield(bvac(1,1), brtemp, extcur(ig), nbvac)
         CALL sum_bfield(bvac(1,2), bptemp, extcur(ig), nbvac)
         CALL sum_bfield(bvac(1,3), bztemp, extcur(ig), nbvac)

      END DO

!
!     MUST ADD EXTERNAL LOOP STUFF LATER
!     MAY DECIDE TO WRITE THAT INFO INTO A SEPARATE FILE
!     FOR NOW, JUST SET DEFAULTS
!
      nobser = 0
      nobd   = 0
      nbsets = 0

      CALL cdf_inquire(ngrid, vn_mgmode, dimlens, ier=istat)
      IF (istat .eq. 0) THEN
         CALL cdf_read(ngrid, vn_mgmode, mgrid_mode)
      ELSE
         mgrid_mode = 'N'
      END IF

      CALL cdf_inquire(ngrid, vn_coilcur, dimlens, ier=istat)
      IF (istat .eq. 0) THEN
	   IF (ALLOCATED(raw_coil_current)) DEALLOCATE(raw_coil_current)
         ALLOCATE (raw_coil_current(nextcur), stat=istat)
         IF (istat .ne. 0) STOP 'Error allocating RAW_COIL in mgrid_mod'
         CALL cdf_read(ngrid, vn_coilcur, raw_coil_current)
      END IF

      CALL cdf_close(ngrid)

      IF (ALLOCATED(brtemp))
     1    DEALLOCATE (vn_br, vn_bz, vn_bp, brtemp, bptemp, bztemp)

      END SUBROUTINE read_mgrid_nc
!DEC$ ENDIF

      SUBROUTINE sum_bfield(bfield, bf_add, cur, n1)
      INTEGER     :: n1
      REAL(rprec), INTENT(inout) :: bfield(n1)
      REAL(rprec), INTENT(in)    :: bf_add(n1)
      REAL(rprec) :: cur

      bfield = bfield + cur*bf_add

      END SUBROUTINE sum_bfield

      SUBROUTINE assign_bptrs(bptr)
      IMPLICIT NONE
      REAL(rprec), TARGET, INTENT(in) :: bptr(nr0b,nz0b,np0b,3)

      brvac => bptr(:,:,:,1)
      bpvac => bptr(:,:,:,2)
      bzvac => bptr(:,:,:,3)

      END SUBROUTINE assign_bptrs
        
      SUBROUTINE free_mgrid (istat)
      INTEGER :: istat 
     
      istat = 0
      mgrid_path_old = ''

      IF (ALLOCATED(bvac)) DEALLOCATE (bvac,stat=istat)
      IF (ALLOCATED(xobser))
     1   DEALLOCATE (xobser, xobsqr, zobser, unpsiext, dsiext,
     2      psiext,plflux, iconnect, needflx, needbfld, plbfld,
     3      nbcoils, rbcoil, zbcoil, abcoil, bcoil, rbcoilsqr, dbcoil,
     4      pfcspec,dsilabel, bloopnames, curlabel, b_chi, stat=istat)
      IF (ALLOCATED(raw_coil_current)) DEALLOCATE(raw_coil_current)

      IF (ALLOCATED(rlim))
     1   DEALLOCATE (rlim,zlim, reslim,seplim,stat=istat)

      END SUBROUTINE free_mgrid

      END MODULE mgrid_mod
