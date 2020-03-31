      MODULE vacfield_mod
      USE stel_kinds
      IMPLICIT NONE
      INTEGER, PARAMETER :: ncur = 100, msmax = 64, nresmax = 10
      INTEGER, PARAMETER :: iout = 10, jout = 11, kout = 12, lout = 13
      INTEGER :: nfp, ngroups, ntmax, npmax, lopt, nres, numres
      INTEGER, DIMENSION(0:nresmax) :: mres, nimax, nrp
      REAL(rprec) :: drsurf, residsum, raxiscon
      REAL(rprec), DIMENSION(msmax) :: ri, zi, di, drho, iota
      REAL(rprec), DIMENSION(0:nresmax) :: vres, rres, dres
      REAL(rprec), DIMENSION(0:nresmax) :: trace, deter, resid
                           !trace, determinant, residue of map
      REAL(rprec), DIMENSION(0:nresmax) :: rres_min, rres_max, reswt
      REAL(rprec), DIMENSION(1:nresmax) :: rmax_isl, rmin_isl
      REAL(rprec), DIMENSION(1:nresmax) :: zmax_isl, zmin_isl
      REAL(rprec), DIMENSION(1:nresmax) :: shear_r, delu_r
      REAL(rprec), DIMENSION(1:nresmax) :: shear_p, delu_p
      REAL(rprec), DIMENSION(1:nresmax) :: area_isl, flux_isl
      REAL(rprec), DIMENSION(1:nresmax) :: alpha
      REAL(rprec), DIMENSION(1:nresmax) :: sigma_resid
      LOGICAL, DIMENSION(1:nresmax) :: lsymm
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: ropnts, zopnts
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: r_isl, z_isl
      LOGICAL, DIMENSION(:,:), ALLOCATABLE :: lli, lri
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: rt, zt
      REAL(rprec), DIMENSION(msmax) :: rmag, zmag
      INTEGER :: ierr, mc, msurf, nvar, nstell_coils
      INTEGER :: ierr_mpi, myid, numprocs
      REAL(rprec), DIMENSION(ncur) :: vac_extcur, xc
      REAL(rprec) :: raxis, zaxis, r0, z0, dsmin, dsurf, tsurf, v0
      REAL(rprec), DIMENSION(nresmax) :: atol, rtol
      LOGICAL :: laxis, lprt, lisize, lsurf, ldisplace, lscreen
      LOGICAL, DIMENSION(:), ALLOCATABLE :: lopnt

      TYPE angles_shifts
         REAL(rprec) :: as_array(3)
      END TYPE angles_shifts


      TYPE (angles_shifts), DIMENSION(20)          :: angles, angles_lo, 
     1     angles_hi, sigma_angles_lo, sigma_angles_hi
      TYPE (angles_shifts), DIMENSION(20)          :: shifts, shifts_lo,
     1     shifts_hi, sigma_shifts_lo, sigma_shifts_hi
      TYPE (angles_shifts), DIMENSION(20,10,10)    :: 
     1     control_pt_init, control_pt_final
      CHARACTER(LEN=100) :: extension, coils_file_extension

      NAMELIST /nlst_vacopt/ vac_extcur, msurf, ntmax, ri, zi, lopt, 
     1   r0, z0, v0, dsmin,  dsurf, tsurf, nfp, laxis, raxis, zaxis, 
     2   vres, rres, rres_min, rres_max, reswt, nrp, dres, mres, numres,
     3   raxiscon, rmin_isl, rmax_isl, zmin_isl, zmax_isl, shear_r,
     4   shear_p, alpha, lsymm, lisize, lsurf, lprt, ldisplace, 
     5   nstell_coils, nvar, rtol, atol, coils_file_extension, 
     6   rmag, zmag, sigma_resid,
     7   angles, angles_lo, angles_hi, sigma_angles_lo, sigma_angles_hi,
     8   shifts, shifts_lo, shifts_hi, sigma_shifts_lo, sigma_shifts_hi,
     9   control_pt_init

      INTERFACE write_invac
      MODULE PROCEDURE write_invac, write_invac2
      END INTERFACE

      CONTAINS

      SUBROUTINE initialize_vacfield (arg1)
      USE safe_open_mod
      USE biotsavart

      INTEGER :: n, np, istat, invac
      CHARACTER*(*), INTENT(in) :: arg1
      CHARACTER :: invac_file*100

!     Get extension of invac file, used for writing output files
      extension = TRIM(arg1)
      n = INDEX(extension, ".", BACK=.TRUE.)
      extension = extension(n+1:)
      coils_file_extension = TRIM(extension)
      
!     Read namelist input file
      invac = 20
      vac_extcur = 0
      invac_file = arg1

      CALL safe_open (invac, istat, invac_file, 'old', 'formatted')
      IF (istat .ne. 0) THEN
         PRINT *, 'Error opening namelist input file', invac_file
         STOP 02
      END IF

      CALL read_vacopt_namelist(invac, istat)

      IF (istat .ne. 0) THEN
         PRINT *,' Error reading ', TRIM(invac_file), 
     1           ' iostat = ', istat
         PRINT *,' Be sure EXTCUR -> VAC_EXTUR!'
         STOP 03
      END IF

      CLOSE(unit=invac)

!     Number of toroidal symmetry planes
      npmax = 2*nfp

!     Read in coils-dot file
      CALL initialize_biotsavart (vac_extcur, coils_file_extension, 
     1                            scaled=.false.)     

      END SUBROUTINE initialize_vacfield

      SUBROUTINE read_vacopt_namelist (iunit, istat)
      INTEGER :: iunit, istat
      INTEGER :: n
!
!     READS IN nls_vacopt NAMELIST
!
      ldisplace = .false.       !.true, will apply rotations, shifts to produce new coils-dot file
      laxis = .true.
      lprt = .false.
      lisize = .false.
      lsurf = .false.
      msurf = 16
      ntmax = 200
      nfp = 2
      lopt = 0
      nvar = 10
      raxis  = 0.7173_dp
      raxiscon = 0.7173_dp
      zaxis  = 0
      rmag(1) = 0.7173_dp
      rmag(2) = 1.2945_dp
      rmag(3) = 0.7173_dp
      rmag(4) = 1.2945_dp
      zmag(1) = 0.0_dp
      zmag(2) = 0.0_dp
      zmag(3) = 0.0_dp
      zmag(4) = 0.0_dp
      r0 = 0.7173_dp
      z0 = 0.0_dp
      v0 = 0.5_dp
      numres = 3                ! number of targeted resonances
      sigma_resid = 1           ! sigma (inverse weight) associated with targeted resonances
      DO n = 1, SIZE(angles,1)
         angles(n) = angles_shifts ( (/0._dp, 0._dp, 0._dp/) )
         shifts(n) = angles_shifts ( (/0._dp, 0._dp, 0._dp/) )
         sigma_angles_lo(n) = angles_shifts ( (/1._dp, 1._dp, 1._dp/) )
         sigma_shifts_lo(n) = angles_shifts ( (/1._dp, 1._dp, 1._dp/) )
         sigma_angles_hi(n) = angles_shifts ( (/1._dp, 1._dp, 1._dp/) )
         sigma_shifts_hi(n) = angles_shifts ( (/1._dp, 1._dp, 1._dp/) )

         !control_pt arrays: fixed pts on coils; init is the initial value (x,y,z)
         !                   final is the rotated+shifted pt (x,y,z) written to chisq file
         control_pt_init(n,:,:) = 
     1                        angles_shifts ( (/0._dp, 0._dp, 0._dp/) )
      END DO
      mres(0) = 1               ! order of fixed point (magnetic axis)
      mres(1) = 10              ! mode number for targeted resonance
      mres(2) = 9               ! mode number for targeted resonance
      mres(3) = 8               ! mode number for targeted resonance
      vres(0) = 0.5_dp          ! toroidal plane for axis search
      vres(1) = 0.0_dp          ! toroidal plane for O-point search
      vres(2) = 0.5_dp          ! toroidal plane for O-point search
      vres(3) = 0.0_dp          ! toroidal plane for O-point search
      rres(0) = 0.72_dp         ! approx. radius of O-point
      rres(1) = 1.33_dp         ! approx. radius of O-point
      rres(2) = 0.90_dp         ! approx. radius of O-point
      rres(3) = 1.37_dp         ! approx. radius of O-point
      nrp(0) = 1                ! no. subintervals in fixed pt. search
      nrp(1) = 1                ! no. subintervals in fixed pt. search
      nrp(2) = 1                ! no. subintervals in fixed pt. search
      nrp(3) = 1                ! no. subintervals in fixed pt. search
      rres_min(0) = 0.70_dp     ! min. R for O-point search
      rres_max(0) = 0.73_dp     ! max. R for O-point search
      rres_min(1) = 1.33_dp     ! min. R for O-point search
      rres_max(1) = 1.36_dp     ! max. R for O-point search
      rres_min(2) = 0.88_dp     ! min. R for O-point search
      rres_max(2) = 0.94_dp     ! max. R for O-point search
      rres_min(3) = 1.37_dp     ! min. R for O-point search
      rres_max(3) = 1.40_dp     ! max. R for O-point search
      reswt(0) = 1.0_dp         ! weight for residue minimization
      reswt(1) = 1.0_dp         ! weight for residue minimization
      reswt(2) = 1.0_dp         ! weight for residue minimization
      reswt(3) = 1.0_dp         ! weight for residue minimization
      dres(0) = 0.010_dp        ! step size for tangent map approx.
      dres(1) = 0.005_dp        ! step size for tangent map approx.
      dres(2) = 0.010_dp        ! step size for tangent map approx.
      dres(3) = 0.005_dp        ! step size for tangent map approx.
      vac_extcur(1) = 0         ! vacuum external current array
      rmin_isl(1) =  1.33_dp    ! min. R for island size calc.
      rmax_isl(1) =  1.36_dp    ! max. R for island size calc.
      zmin_isl(1) = -0.40_dp    ! min. Z for island size calc.
      zmax_isl(1) =  0.40_dp    ! max. Z for island size calc.
      rmin_isl(2) =  0.88_dp    ! min. R for island size calc.
      rmax_isl(2) =  0.94_dp    ! max. R for island size calc.
      zmin_isl(2) = -0.10_dp    ! min. Z for island size calc.
      zmax_isl(2) =  0.10_dp    ! max. Z for island size calc.
      rmin_isl(3) =  1.37_dp    ! min. R for island size calc.
      rmax_isl(3) =  1.40_dp    ! max. R for island size calc.
      zmin_isl(3) = -0.40_dp    ! min. Z for island size calc.
      zmax_isl(3) =  0.40_dp    ! max. Z for island size calc.
      shear_r(1) =   0.56_dp    ! shear (1/m) of unperturbed flow at res.
      shear_r(2) =   0.56_dp    ! shear (1/m) of unperturbed flow at res.
      shear_r(3) =   0.53_dp    ! shear (1/m) of unperturbed flow at res.
      shear_p(1) =  -0.3875_dp  ! shear (1/Wb) of unperturbed flow at res.
      shear_p(2) =  -0.3875_dp  ! shear (1/Wb) of unperturbed flow at res.
      shear_p(3) =  -0.2813_dp  ! shear (1/Wb) of unperturbed flow at res.
      lsymm(1) = .false.        ! = true if starting in elongated plane
      lsymm(2) = .false.        ! = true if starting in elongated plane
      lsymm(3) = .false.        ! = true if starting in elongated plane
      alpha(1) = 1.0_dp         ! factor to adjust du(m) for starting pt.
      alpha(2) = 1.0_dp         ! factor to adjust du(m) for starting pt.
      alpha(3) = 1.0_dp         ! factor to adjust du(m) for starting pt.
      dsmin = 0.00_dp
      dsurf = 0.01_dp
      tsurf = 0.00_dp
      rtol = 1.e-8_dp
      atol = 1.e-8_dp

      DO n = 1, msmax
         di(n) = dsmin + n*dsurf
         ri(n) = raxis + di(n)*COS(tsurf)
         zi(n) = zaxis + di(n)*SIN(tsurf)
      END DO

      nstell_coils = -1       !No. of unique "stellarator-symmetric" coils to rotate
      
      READ (iunit, nml=nlst_vacopt, iostat=istat)

      END SUBROUTINE read_vacopt_namelist

      SUBROUTINE write_invac (iunit, istat)
      USE stel_constants
      USE write_array_generic
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: iunit, istat
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, j, k
      CHARACTER(LEN=100) :: form1, form2
!-----------------------------------------------

      WRITE(iunit,'(a12)')'&NLST_VACOPT'
      WRITE(iunit,100) 'LDISPLACE = ', ldisplace
      WRITE(iunit,100) 'LAXIS = ', laxis
      WRITE(iunit,100) 'LPRT = ', lprt
      WRITE(iunit,100) 'LISIZE = ', lisize
      WRITE(iunit,200) 'LOPT = ', lopt
      WRITE(iunit,200) 'MSURF = ', msurf
      WRITE(iunit,200) 'NTMAX = ', ntmax
      WRITE(iunit,200) 'NFP = ', nfp
      WRITE(iunit,200) 'NSTELL_COILS = ', nstell_coils
      WRITE(iunit,200) 'NVAR = ', nvar
      WRITE(iunit,300) 'RAXIS = ', raxis
      WRITE(iunit,300) 'R0 = ', r0
      WRITE(iunit,300) 'Z0 = ', z0
      WRITE(iunit,300) 'V0 = ', v0
      WRITE(iunit,300) 'DSURF = ', dsurf
      WRITE(iunit,300) 'DSMIN = ', dsmin
      WRITE(iunit,300) 'TSURF = ', tsurf
      WRITE(iunit,200) 'NUMRES = ', numres
      WRITE(iunit,'(1x,3a)') "COILS_FILE_EXTENSION = '",
     1                        TRIM(coils_file_extension),"'"
      IF (msmax .gt. 0) THEN
         CALL write_array(iunit,'RMAG', rmag, msmax)
         CALL write_array(iunit,'ZMAG', zmag, msmax)
         CALL write_array(iunit,'RI', ri, msurf)
         CALL write_array(iunit,'ZI', zi, msurf)
      END IF
      IF (nresmax .gt. 0) THEN
         CALL write_array(iunit,'MRES', mres, nresmax+1, low_index=0)
         CALL write_array(iunit,'VRES', vres, nresmax+1, low_index=0)
         CALL write_array(iunit,'RRES', rres, nresmax+1, low_index=0)
         CALL write_array(iunit,'NRP', nrp, nresmax+1, low_index=0)
         CALL write_array(iunit,'RRES_MIN', rres_min, nresmax+1, 
     1                    low_index=0)
         CALL write_array(iunit,'RRES_MAX', rres_max, nresmax+1,
     1                    low_index=0)
         CALL write_array(iunit,'RESWT', reswt, nresmax+1, low_index=0)
         CALL write_array(iunit,'DRES', dres, nresmax+1, low_index=0)
      END IF
      IF (ncur .gt. 0) THEN
         CALL write_array(iunit,'VAC_EXTCUR', vac_extcur, ncur)
      END IF
      IF (nresmax .gt. 0) THEN
         CALL write_array(iunit,'RMIN_ISL', rmin_isl, nresmax)
         CALL write_array(iunit,'RMAX_ISL', rmax_isl, nresmax)
         CALL write_array(iunit,'ZMIN_ISL', zmin_isl, nresmax)
         CALL write_array(iunit,'ZMAX_ISL', zmax_isl, nresmax)
         CALL write_array(iunit,'SHEAR_R', shear_r, nresmax)
         CALL write_array(iunit,'SHEAR_P', shear_p, nresmax)
         CALL write_array(iunit,'LSYMM', lsymm, nresmax)
         CALL write_array(iunit,'ALPHA', alpha, nresmax)
         CALL write_array(iunit,'SIGMA_RESID', sigma_resid, nresmax)
      END IF

      DO i = 1, nstell_coils
         IF (i .lt. 10) THEN
            form1 = '(2x,a,i1,a,1p3e15.6)'
            form2 = '(2x,a,i1,2(a1,i2.2),a,1p3e15.6)'
         ELSE 
            form1 = '(2x,a,i2,a,1p3e15.6)'
            form2 = '(2x,a,i2,2(a1,i2.2),a,1p3e15.6)'
         END IF
         WRITE(iunit, form1, iostat=istat)
     1       'ANGLES(',i,') = ', angles(i)
         WRITE(iunit, form1, iostat=istat)
     1       'SHIFTS(',i,') = ', shifts(i)
         WRITE(iunit, form1, iostat=istat)
     1       'ANGLES_LO(',i,') = ', angles_lo(i)
         WRITE(iunit, form1, iostat=istat)
     1       'SHIFTS_LO(',i,') = ', shifts_lo(i)
         WRITE(iunit, form1, iostat=istat)
     1       'ANGLES_HI(',i,') = ', angles_hi(i)
         WRITE(iunit, form1, iostat=istat)
     1       'SHIFTS_HI(',i,') = ', shifts_hi(i)
         WRITE(iunit, form1, iostat=istat)
     1       'SIGMA_ANGLES_LO(',i,') = ', sigma_angles_lo(i)
         WRITE(iunit, form1, iostat=istat)
     1       'SIGMA_SHIFTS_LO(',i,') = ', sigma_shifts_lo(i)
         WRITE(iunit, form1, iostat=istat)
     1       'SIGMA_ANGLES_HI(',i,') = ', sigma_angles_hi(i)
         WRITE(iunit, form1, iostat=istat)
     1       'SIGMA_SHIFTS_HI(',i,') = ', sigma_shifts_hi(i)

         DO j = 1, SIZE(control_pt_init,2)
            DO k = 1, SIZE(control_pt_init,3)
               WRITE(iunit, form2, iostat=istat)
     1            'CONTROL_PT_INIT(',i,',',j,',',k,') = ', 
     2             control_pt_init(i,j,k)
            END DO
         END DO

      END DO

      WRITE(iunit,'(a)')'/'

 100  FORMAT(4(1x,a,l2,','))
 200  FORMAT(4(1x,a,i6,','))
 300  FORMAT(3(1x,a,1pe12.4,','))
 350  FORMAT(2(1x,a,1pe12.4,','))
 450  FORMAT(1x,a,1pe12.4,',')

      END SUBROUTINE write_invac

      SUBROUTINE write_invac2 (filename, istat)
      USE stel_constants
      USE write_array_generic
      USE safe_open_mod
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: iunit, istat
      CHARACTER(LEN=*) :: filename
!-----------------------------------------------
      CALL safe_open(iunit, istat, filename, 'replace', 'formatted')
      IF (istat .ne. 0) RETURN

      CALL write_invac (iunit, istat)

      CLOSE (iunit)

      END SUBROUTINE write_invac2

      END MODULE vacfield_mod
