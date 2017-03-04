      MODULE PERTURBATION
      USE quantities
      USE descriptor_mod, ONLY: iam
      USE siesta_state, ONLY: update_state, update_state_par
      
      INTEGER  :: nsin=101, mpolin=12, ntorin=3
      INTEGER  :: mres1, nres1, mres2, nres2                          !Backwards compat
      INTEGER  :: mres(10)
      INTEGER  :: niter_max
      REAL(dp) :: HelPert1, HelPert2, rad1, rad2                      !Backwards compat
      REAL(dp) :: HelPert(10)
      REAL(dp) :: ftol = 1.E-20_dp
      REAL(dp) :: eta_factor = 0.1_dp
      LOGICAL  :: ladd_pert = .TRUE.
      LOGICAL  :: lresistive = .TRUE.                                 ! RS: If true, add resistive electric field.
      LOGICAL  :: lrestart = .FALSE.
      LOGICAL  :: lprecon_diag = .TRUE.
      LOGICAL  :: lposDef=.FALSE.                          !<use pos-def approx Hessian if TRUE for fsq > 10^-12
      CHARACTER(LEN=256) :: wout_file, restart_ext
      REAL(dp), ALLOCATABLE :: buv_res(:,:,:)

!NETCDF PARAMETERS
      CHARACTER(LEN=*), PARAMETER :: vn_nsin='nrad', vn_mpolin='mpol',  &
         vn_ntorin='ntor', vn_wout='wout file', vn_bsups='Bsups(m,n,r)',&
         vn_bsupu='Bsupu(m,n,r)', vn_bsupv='Bsupv(m,n,r)',              &
         vn_pres='pres(m,n,r)'
      CHARACTER(LEN=*), DIMENSION(3), PARAMETER ::                      &
                  r3dim = (/'m-mode','n-mode','radius'/)

!SCRATCH ARRAYS
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: temp1, temp2, temp3, temp4

!OUTPUT FILE PARAMETERS
      INTEGER:: nrs=200, nzs=200, nphis=2                             ! BField resolution (CYL):  R and Z
      INTEGER:: nvs=100, nus=150, nss=150                             ! BField resolution (VMEC): V (and PHI), U, S      
      LOGICAL :: l_tracing = .FALSE.                                  ! Do not produce file with CYL field for Poincare plots.
      LOGICAL :: l_silo_output = .FALSE.
      LOGICAL :: l_silo3D = .FALSE.
      LOGICAL :: l_output_alliter = .FALSE.                           ! Default: do not dump output at every iteration on file (RS, Jul'10)
      LOGICAL :: l_VMEC_Uniform

      CHARACTER(LEN=256) :: filename

      CONTAINS

      SUBROUTINE init_data (nprecon)
      USE date_and_computer
      USE Hessian, ONLY: levmarq_param, mupar, levmarq_param0, mupar0
      IMPLICIT NONE
      INTEGER, INTENT(out) :: nprecon
      INTEGER              :: istat, index1, imon, niter
      CHARACTER(LEN=10)    :: date0, time0, zone0
      CHARACTER(LEN=256)   :: temp, short_name
      NAMELIST /SIESTA_INFO/ ladd_pert, lresistive, lrestart, niter,   &
               nsin, mpolin, ntorin, ftol, wout_file, restart_ext,     &
               mres1, nres1, rad1, helpert1, mres2, nres2, rad2, helpert2, &
               mres, HelPert, levmarq_param, mupar, eta_factor,        &
               nphis, nrs, nzs, nvs, nus, nss, l_tracing,              &
               lprecon_diag, lPosDef,                                  &
               l_silo_output, l_silo3D, nprecon, l_output_alliter      & ! If l_output_alliter, dump info (SILO or TRACING) on file at every iteration (RS, Jul'10)
               ,l_VMEC_Uniform

!-----------------------------------------------
!     
!     READS siesta.jcf (JOB CONTROL FILE) for data
!
      levmarq_param = 1.E-3_dp
      mupar         = 0
      niter         = 10                !Number iterations after diagonal prec
      HelPert1      = 0
      HelPert2      = 0
      mres          = 0
      HelPert       = 0
      lresistive    = .TRUE.

      OPEN (unit=33, file='siesta.jcf', delim='apostrophe', iostat=istat)
      IF (istat .ne. 0) STOP 'SIESTA.JCF file not found'

      READ (33, NML=SIESTA_INFO, iostat=istat) 
      IF (istat .ne. 0) THEN
         IF (iam .EQ. 0) THEN
            STOP 'Error reading siesta.jcf file'
         ELSE
            STOP
         END IF
      END IF
      CLOSE (unit=33)

      nss = nsin
      levmarq_param0 = levmarq_param
      mupar0         = mupar
      niter_max      = niter

      temp = wout_file
      index1 = INDEX(temp,'WOUT')
      IF (index1 .eq. 0) index1 = INDEX(temp,'wout')
      IF (index1 .eq. 0) STOP 'ERROR: WRONG WOUT FILE NAME'
      temp = wout_file(index1:)
!
!     STRIP FILE EXTENSION: NAME COULD INCLUDE '.' IN IT!
!
      index1 = SCAN(temp, '.',BACK=.TRUE.)
      IF (index1 == 0) THEN
         temp = TRIM(temp(6:))
      ELSE
         temp = TRIM(temp(6:index1-1))
      END IF

      IF (HelPert1.ne.zero .or. HelPert2.ne.zero) THEN
         HelPert(1) = HelPert1; mres(1) = mres1
         HelPert(2) = HelPert2; mres(2) = mres2
      END IF

      IF (iam .ne. 0) RETURN

!     WOUT FILE MUST BE OF FORM "wout_NAME.nc (or txt)" or "wout_NAME"
      WRITE (filename,'(3a,i4.4,a,2(i3.3,a))') "output_", TRIM(temp),   &
             '_',nsin,'X',mpolin,'X',ntorin,".txt"
      OPEN (unit=33, file=filename,iostat=istat)
      IF (istat .ne. 0) STOP 'ERROR WRITING SIESTA OUTPUT FILE'

      DO istat = 6,33,27
	     WRITE (istat, 10)                                          & 
               'SIESTA MHD EQUILIBRIUM CODE v2.3 (050213)',             &
               'Scalable Island Equilibrium ',                          &
               'Solver for Toroidal Applications' 
         IF (HelPert1.ne.zero .or. HelPert2.ne.zero) THEN
            WRITE (istat, *) 'Please update jcf file to new format'
            WRITE (istat, *) 'Replace mres1,2 with mres array, HelPert1,2 with HelPert array'
         END IF
      END DO

10    FORMAT (62('-'),/,1x,a,/,2(1x,a),/,62('-'),/)

      WRITE (33, *)'CASE: ', TRIM(temp)
      PRINT *,'OUTPUT FILE (SCREEN DUMP): ', "output_" // TRIM(temp) // ".txt"

      CALL DATE_AND_TIME(date0,time0,zone0)
      READ (date0(5:6),'(i2)')imon
      WRITE (33,20) months(imon),date0(7:8),date0(1:4),                 &
     &              time0(1:2),time0(3:4),time0(5:6)

20    FORMAT(' DATE = ',a3,' ',a2,',',a4,' ',' TIME = ',2(a2,':'),a2,/)

      DO istat = 1, SIZE(HelPert)
         IF (HelPert(istat) .ne. zero) THEN
            WRITE (33, 30) istat, mres(istat), HelPert(istat)
         END IF
      END DO

30    FORMAT(i2,' mres: ',i4,' HelPert: ',1pe9.2)

      END SUBROUTINE init_data

      SUBROUTINE add_perturb (xc, getwmhd)
      USE nscalingtools, ONLY: startglobrow, endglobrow
      IMPLICIT NONE
      INTEGER     :: neqs
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp)            :: xc(:)
      REAL(dp), EXTERNAL  :: getwmhd
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER,PARAMETER  :: irad_scale=10, ipert_sign=2
      INTEGER            :: js, l, iprint, mres0, nres0, isign, icount, &
                            irscale, imin(2), jstart, istat, jsave
      INTEGER            :: nsmin, nsmax
      REAL(rprec)        :: wmhd_array(irad_scale,ipert_sign), eta_factor_save
      REAL(rprec)        :: w0, wmhd, normal, HelPert0, rad, chip0, phip0
      CHARACTER(LEN=1)   :: chsign
      LOGICAL            :: lresist_save, lpar
!-----------------------------------------------
!
!     Add helical RESISTIVE flux perturbation on full mesh
!     Approximate nonideal E|| ~ f(mu+nv) B with B ~ Bsubv
!
      IF (ntor .eq. 0) RETURN

!     INITIAL STATE
 
      lresist_save = lresistive
      lresistive=.FALSE.
      eta_factor_save = eta_factor
      eta_factor=1
      xc = 0
#if defined(SKS)
      lpar=.TRUE.
      nsmin=MAX(1,startglobrow-1); nsmax=MIN(endglobrow+1,ns)
#else
      lpar = .FALSE.
#endif
      w0 = getwmhd(xc, lpar)
      wmhd_array = 2*w0                                                   !DON'T WANT TO COUNT THIS UNLESS IT IS ACTUALLY A RESONANT MODE

      ALLOCATE (buv_res(ntheta,nzeta,ns), stat=istat)
      IF (istat .NE. 0) STOP 'ALLOCATION ERROR IN add_perturbation'

      IF (iam .eq. 0) THEN
      DO iprint = 6,33,27
         WRITE (iprint, '(/,a,/,2a)')                                    &
         ' Adding helical flux perturbations',                           &
         '  Phase   10^6 X Del-W    mres    nres    HelPert',            &
         '     rad   |m*chip+n*phip|'
      END DO
      ENDIF

!     Resonant resistitive perturbations first
      lresistive=.TRUE.
      normal = hs_i**2

      RES_TEST: DO icount = 1, SIZE(HelPert)
         HelPert0 = ABS(HelPert(icount)/normal)
         mres0 = mres(icount)
         IF (HelPert0.EQ.zero .or. mres0.GT.mpol) CYCLE

!Scan in radius to determine primary resonance
         jstart = 4

 100     CONTINUE
         CALL FindResonance(mres0, nres0, rad, jstart)
!SPH: OFF 030314         IF (rad .EQ. zero) CYCLE
         IF (rad .EQ. zero) CYCLE

      PERT_SIGN: DO isign = 1, 2
#if defined(SKS)
         lpar=.TRUE.
#endif      

!     NOTE: in update_bfield, where perturbed b-field is computed, this is multiplied
!           by eta_prof = rho*(1-rho) so it vanishes at both ends. Multiply again
!           here by that factor to assure radial derivatives entering in dB also vanish
!           at endpoints r=0,1
         SCALE_TEST: DO irscale = 1, irad_scale
            CALL GetResPert(irscale, isign, mres0, nres0, rad, HelPert0, chip0, phip0)
            xc = 0
            wmhd = getwmhd(xc, lpar)
            wmhd_array(irscale,isign) = wmhd

!            IF (iam.eq.0) PRINT *,'isign: ',isign, ' irscale: ',irscale,' WMHD: ',wmhd

         END DO SCALE_TEST
      END DO PERT_SIGN

!     Make sure energy decreases
      imin = MINLOC(wmhd_array)
      wmhd = wmhd_array(imin(1),imin(2))
!      IF (wmhd < w0) THEN
!     recompute perturbation buv_rex here
         irscale = imin(1)
         isign   = imin(2)
         IF (iam .eq. 0) THEN
         IF (isign == 2) THEN
            chsign = '-'
         ELSE
            chsign = '+'
         END IF

         DO iprint = 6, 33, 27
            WRITE (iprint, '(4x,a,1p,e16.3,2(i8),3x,1pe9.2,0p,f8.2,f11.2)') &
                chsign,10**6*(wmhd-w0), mres0, nres0, HelPert0*normal,  &
                rad, ABS(mres0*chip0+nres0*nfp*phip0)
         END DO
         END IF
         CALL GetResPert(irscale, isign, mres0, nres0,                  &
                         rad, HelPert0, chip0, phip0)
         xc = 0
!         lpar = .FALSE.
         wmhd = getwmhd(xc, lpar)
         IF (iam.EQ.0 .AND. ABS(wmhd-wmhd_array(irscale,isign)).GT. 1.E-12_dp) STOP 'Error1 in Perturb'
         IF (lpar) THEN
            CALL update_state_par(.FALSE., zero, zero)
         ELSE
            CALL update_state(.FALSE., zero, zero)
         END IF
         w0 = wmhd
!      END IF
 
 1001 CONTINUE

         IF (jstart .lt. (ns-10)) THEN
            jstart = jstart+10             !!Avoid flat iota multiples
            GOTO 100                       !!Look for multiple resonances
         END IF

      END DO RES_TEST

      IF (iam .eq. 0) THEN
         DO iprint = 6, 33, 27
            WRITE (iprint, *)
         END DO
      END IF

      DEALLOCATE (buv_res, stat=istat)
      lresistive = lresist_save
      eta_factor = eta_factor_save

      END SUBROUTINE add_perturb


      SUBROUTINE FindResonance(mres0, nres0, resrad, jstart)
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in)      :: mres0
      INTEGER, INTENT(out)     :: nres0
      INTEGER, INTENT(inout)   :: jstart
      REAL(rprec), INTENT(out) :: resrad
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER                   :: js, nmin, jmin, nres02
      REAL(rprec)               :: del1, del2, rho1, rho2, delmin
      LOGICAL                   :: lZero
!-----------------------------------------------
      del2 = 0
      resrad = 0
      delmin = -1
      lZero = .FALSE.
      DO js = jstart, ns-1
         del1 = mres0*(chipf_i(js)/nfp)
         IF (del1 .eq. zero) CYCLE
         IF (phipf_i(js) .ne. zero) THEN
            nres0 = -NINT(del1/phipf_i(js))
         ELSE 
            nres0 = 0
         ENDIF
         del1 = del1 + nres0*phipf_i(js)
         del1 = del1/MAX(ABS(chipf_i(js)),ABS(phipf_i(js)),1.E-3_dp)
         IF (delmin.eq.-1 .or. ABS(del1).lt.delmin) THEN
            delmin = ABS(del1)
            nmin = nres0
            jmin = js
         END IF
         
         IF (del2 .ne. zero) THEN     !!adjust del2 to SAME value nres0
            del2 = del2 + (nres0-nres02)
         END IF
            
         IF (del1*del2 < zero) THEN
            lZero = .TRUE.
            jstart = js+1
            rho2 = hs_i*(js-2)
            rho1 = hs_i*(js-1)
            resrad = (rho2*ABS(del1)+rho1*ABS(del2))/(ABS(del1)+ABS(del2))
            delmin=0
            EXIT
         END IF
         del2 = del1
         nres02 = nres0
      END DO 

      IF (delmin.lt.ABS(nmin)*1.E-2_dp .AND. ABS(nmin).le.ntor) THEN
         jstart = jmin
         nres0 = nmin
         resrad = hs_i*(jmin-1)
      ELSE
         jstart=ns-1
         resrad = 0
      END IF

!Do not search for multiple roots if no crossings occurred
      IF (.not.lZero) jstart=ns-1

      END SUBROUTINE FindResonance


      SUBROUTINE GetResPert(irscale, isign, mres0, nres0, rad,          &
                            HelPert0, chip0, phip0)
      USE fourier, ONLY: toijsp2
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in)       :: irscale, isign
      INTEGER, INTENT(in)       :: mres0, nres0
      REAL(rprec), INTENT(in)   :: rad
      REAL(rprec), INTENT(out)  :: chip0, phip0
      REAL(rprec), INTENT(in)   :: HelPert0
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER     :: js, istat
      REAL(rprec) :: rho, rhores, pert_prof(ns), HelP_local
      REAL(rprec), ALLOCATABLE   :: bmn_res(:,:)
!-----------------------------------------------
      IF (mres0.gt.mpolin .or. ABS(nres0).gt.ntorin) RETURN

      HelP_local = HelPert0
      IF (isign == 2) HelP_local = -HelPert0

      ALLOCATE (bmn_res(0:mpol,-ntor:ntor), stat=istat)
      IF (istat .ne. 0) STOP 'ISTAT != 0 IN GETRESPERT'
      bmn_res = 0

      js = INT(rad/hs_i) + 1
      IF (js.LT.1 .OR. js.GT.ns) RETURN
      chip0 = chipf_i(js)
      phip0 = phipf_i(js)
      DO js = 1,ns
         rho = hs_i*(js-1)
         rhores = rho-rad
         rhores = irscale*rhores                      !decay length: about 1/irscale of minor radius
         pert_prof(js) = one/(one+rhores*rhores)      !tearing parity
!         pert_prof(js) = rhores*pert_prof(js)        !odd-parity (rippling mode)
      END DO

      pert_prof = HelP_local*pert_prof
         
      bmn_res(mres0,nres0) = 1
      CALL toijsp2(bmn_res, buv_res, 0, 1)
      DO js = 2,ns
         buv_res(:,:,js) = buv_res(:,:,1)*pert_prof(js)
      END DO

      buv_res(:,:,1) = pert_prof(1)*buv_res(:,:,1)

      DEALLOCATE (bmn_res, stat=istat)

      END SUBROUTINE GetResPert
        

      SUBROUTINE write_restart_file
      USE quantities, ONLY: jbsupsmnsh, jbsupumnch, jbsupvmnch, jpmnch, &
                            ns, mpol, ntor
      USE descriptor_mod, ONLY: iam
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER           :: istat, js, m, n
      CHARACTER*(256)   :: filename
!-----------------------------------------------
    
      IF (iam .ne. 0) RETURN

#if defined(NETCDF)
      CALL write_restart_ncfile (istat)
      IF (istat .EQ. 0) RETURN
#endif

      filename = "siesta_" // TRIM(restart_ext) // ".txt"
      OPEN (unit=15, file=TRIM(filename), status='replace',             &
            form='formatted', iostat=istat)
      IF (istat .NE. 0) STOP 'ERROR WRITING SIESTA RESTART FILE'

!     PRINT *,'WRITING RESTART FILE: ', TRIM(filename)
!
!     USE THIS DIMENSION INFO TO COMPARE WITH MESH IN JOB-CONTROL FILE
!     EVENTUALLY USE FOR MULTIGRIDDING
!
      WRITE (15,100) nsin, mpolin, ntorin
      WRITE (15,200) wout_file

      DO js=1,nsin
      DO m =0,mpolin
      DO n = -ntorin,ntorin
      WRITE (15,300) jbsupsmnsh(m,n,js), jbsupumnch(m,n,js),            &
     &               jbsupvmnch(m,n,js), jpmnch(m,n,js)
      ENDDO
      ENDDO
      ENDDO

      CLOSE (unit=15)

 100  FORMAT(3i8)
 200  FORMAT(a256)
 300  FORMAT(1p4e24.16)

      END SUBROUTINE write_restart_file
     
      
#if defined(NETCDF)
      SUBROUTINE write_restart_ncfile (istat)
      USE quantities, ONLY: jbsupsmnsh, jbsupumnch, jbsupvmnch, jpmnch, &
                            ns, mpol, ntor
      USE ezcdf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: istat
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER           :: ncdf, js, m, n
      CHARACTER*(256)   :: filename
!-----------------------------------------------
      filename = "siesta_" // TRIM(restart_ext) // ".nc"

      CALL cdf_open(ncdf,filename,'w',istat)
      IF (istat .NE. 0) RETURN

!     PRINT *,'WRITING RESTART FILE: ', TRIM(filename)

!================================
! Define Variables
!================================
!  Scalars 
      CALL cdf_define(ncdf, vn_nsin, nsin)
      CALL cdf_define(ncdf, vn_mpolin, mpolin)
      CALL cdf_define(ncdf, vn_ntorin, ntorin)
      CALL cdf_define(ncdf, vn_wout, wout_file)

!  3D arrays
      CALL cdf_define(ncdf, vn_bsups, jbsupsmnsh, dimname=r3dim)
      CALL cdf_define(ncdf, vn_bsupu, jbsupumnch, dimname=r3dim)
      CALL cdf_define(ncdf, vn_bsupv, jbsupvmnch, dimname=r3dim)
      CALL cdf_define(ncdf, vn_pres,  jpmnch,     dimname=r3dim)

!================================
! Write Variables
!================================
      CALL cdf_write(ncdf, vn_nsin, nsin)
      CALL cdf_write(ncdf, vn_mpolin, mpolin)
      CALL cdf_write(ncdf, vn_ntorin, ntorin)
      CALL cdf_write(ncdf, vn_wout, wout_file)

      CALL cdf_write(ncdf, vn_bsups, jbsupsmnsh)
      CALL cdf_write(ncdf, vn_bsupu, jbsupumnch)
      CALL cdf_write(ncdf, vn_bsupv, jbsupvmnch)
      CALL cdf_write(ncdf, vn_pres,  jpmnch)

      CALL cdf_close(ncdf)

      END SUBROUTINE write_restart_ncfile
#endif

      SUBROUTINE read_restart_file (nprecon)
      USE quantities, ONLY: jbsupsmnsh, jbsupumnch, jbsupvmnch, jpmnch
      IMPLICIT NONE
      INTEGER, INTENT(inout) :: nprecon
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: m0=0, m1=1
      INTEGER            :: istat, m, n, js, ns, ntor, mpol, nmin
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: temp
      REAL(rprec)        :: r0
      CHARACTER*(256)    :: filename
!-----------------------------------------------
#if defined(NETCDF)
      CALL read_restart_ncfile (istat, ns, ntor, mpol)
      IF (istat .EQ. 0) GOTO 111
#endif
      filename = "siesta_" // TRIM(restart_ext) // ".txt"
      OPEN (unit=15, file=TRIM(filename), status='old',                 &
            form='formatted', iostat=istat)
      IF (istat .NE. 0) THEN 
         PRINT *,'SIESTA RESTART FILE WAS NOT READ: STARTING AT BEGINNING'
         nprecon = 0
         RETURN
      END IF

      IF (iam .eq. 0) THEN
         PRINT *,'READING RESTART FILE: ', TRIM(filename)
      ENDIF
!
!     USE THIS DIMENSION INFO TO COMPARE WITH MESH IN JOB-CONTROL FILE
!     EVENTUALLY USE FOR MULTIGRIDDING
!
      READ (15,100) ns, mpol, ntor
      ALLOCATE (temp1(0:mpol,-ntor:ntor,ns), temp2(0:mpol,-ntor:ntor,ns), &
                temp3(0:mpol,-ntor:ntor,ns), temp4(0:mpol,-ntor:ntor,ns), &
                stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error in read_restart_file'

      READ (15,200) wout_file

      DO js=1,ns
      DO m=0,mpol
      DO n=-ntor,ntor
      READ (15,300) temp1(m,n,js), temp2(m,n,js),                       &
                    temp3(m,n,js), temp4(m,n,js)
      ENDDO
      ENDDO
      ENDDO

 111  CONTINUE

      nmin = MIN(ntor, ntorin)
      jbsupsmnsh(:,:,1) = 0
      jbsupvmnch(:,:,1) = 0
      jpmnch(:,:,1) = 0

      CALL interpit(temp1, jbsupsmnsh, ns, nsin, mpol, mpolin, ntor, ntorin)
      jbsupsmnsh(m1,-nmin:nmin,1) = jbsupsmnsh(m1,-nmin:nmin,2)
      CALL interpit(temp2, jbsupumnch, ns, nsin, mpol, mpolin, ntor, ntorin)
      jbsupumnch(m1,-nmin:nmin,1) = jbsupumnch(m1,-nmin:nmin,2)
      CALL interpit(temp3, jbsupvmnch, ns, nsin, mpol, mpolin, ntor, ntorin)
      jbsupvmnch(m0,-nmin:nmin,1) = jbsupvmnch(m0,-nmin:nmin,2)
      CALL interpit(temp4, jpmnch,     ns, nsin, mpol, mpolin, ntor, ntorin)
      jpmnch(m0,-nmin:nmin,1) = jpmnch(m0,-nmin:nmin,2)

      DEALLOCATE (temp1, temp2, temp3, temp4, stat=istat)
      CLOSE (unit=15)

 100  FORMAT(3i8)
 200  FORMAT(a256)
 300  FORMAT(1p4e24.16)

      END SUBROUTINE read_restart_file

#if defined(NETCDF)
      SUBROUTINE read_restart_ncfile (istat, ns, ntor, mpol)
      USE ezcdf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: istat, ns, mpol, ntor
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER           :: ncdf
      CHARACTER*(256)   :: filename
!-----------------------------------------------
      filename = "siesta_" // TRIM(restart_ext) // ".nc"

      CALL cdf_open(ncdf,filename,'r',istat)
      IF (istat .ne. 0) RETURN

      CALL cdf_read(ncdf, vn_nsin, ns)
      CALL cdf_read(ncdf, vn_mpolin, mpol)
      CALL cdf_read(ncdf, vn_ntorin, ntor)
      CALL cdf_read(ncdf, vn_wout, wout_file)

      ALLOCATE (temp1(0:mpol,-ntor:ntor,ns), temp2(0:mpol,-ntor:ntor,ns), &
                temp3(0:mpol,-ntor:ntor,ns), temp4(0:mpol,-ntor:ntor,ns), &
                stat=istat)
      IF (istat .ne. 0) STOP 'Memory allocation error in read_restart_ncfile'

!================================
! Read Variables
!================================

      CALL cdf_read(ncdf, vn_bsups, temp1)
      CALL cdf_read(ncdf, vn_bsupu, temp2)
      CALL cdf_read(ncdf, vn_bsupv, temp3)
      CALL cdf_read(ncdf, vn_pres,  temp4)

      CALL cdf_close(ncdf)

      END SUBROUTINE read_restart_ncfile
#endif
      
      SUBROUTINE interpit(aold, anew, ns_old, ns_new, m1, m2, n1, n2)
      USE stel_constants, ONLY: one
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in)      :: n1, n2, m1, m2, ns_old, ns_new
      REAL(rprec), INTENT(in)  :: aold(0:m1,-n1:n1,ns_old)
      REAL(rprec), INTENT(out) :: anew(0:m2,-n2:n2,ns_new)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), PARAMETER   :: p5 = 0.5_dp
      REAL(rprec), ALLOCATABLE :: temp1(:,:,:)
      INTEGER                  :: m12, n12, js_old, js_new, istat
      REAL(rprec)              :: hs_old, hs_new, snew, sold, x1
!-----------------------------------------------

      m12 = MIN(m1,m2);  n12 = MIN(n1, n2)
      
      anew = 0

      IF (ns_old .EQ. ns_new) THEN
         anew(0:m12,-n12:n12,:) = aold(0:m12,-n12:n12,:)
         RETURN
      END IF

      ALLOCATE (temp1(0:m12,-n12:n12,ns_old), stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error #1 in INTERPIT'
!
!     FIRST INTERP TO FULL MESH FROM HALF-MESH
!     ASSUMES js=1 FULL-MESH VALUE VALUE IS STORED IN aold(1)
!
      temp1(:,:,1) = 0
      temp1(:,:,2:ns_old-1) = p5*(aold(:,:,2:ns_old-1) + aold(:,:,3:ns_old))
      temp1(:,:,ns_old) = 2*aold(:,:,ns_old) - temp1(:,:,ns_old-1)
      temp1(0:1,:,1)  = 2*aold(0:1,:,2) - temp1(0:1,:,2)
      hs_old = one/(ns_old-1)
      hs_new = one/(ns_new-1)

      DO js_new=1, ns_new
         snew = (js_new-1)*hs_new
         js_old = 1 + INT(snew/hs_old)
         x1 = snew - (js_old-1)*hs_old
         x1 = x1/hs_old
         IF (x1 .GT. one) x1 = 1
         IF (x1 .LT. zero) x1 = 0
         IF (js_old .gt. ns_old) STOP 'js_old >= ns_old in INTERPIT'
         IF (js_old .eq. ns_old) THEN
             anew(:,:,js_new) = temp1(:,:,js_old)
         ELSE
             anew(:,:,js_new) = (one-x1)*temp1(:,:,js_old)     &
                              + x1*temp1(:,:,js_old+1)
         END IF
      END DO
!      
!     AVERAGE BACK TO HALF-MESH
!
      DEALLOCATE (temp1)
      ALLOCATE (temp1(0:m12,-n12:n12,ns_new), stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error #2 in INTERPIT'
      temp1 = anew
      anew(:,:,1) = 0
      anew(:,:,2:ns_new) = p5*(temp1(:,:,2:ns_new) + temp1(:,:,1:ns_new-1))
      DEALLOCATE(temp1)

      END SUBROUTINE interpit

      SUBROUTINE Compute_Resonant_Epar(xc, epar)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(out) :: epar(ns,ntheta,nzeta)
      REAL(dp), INTENT(in)  :: xc(ns,0:mpol,-ntor:ntor,3)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), ALLOCATABLE :: eparmn(:,:,:)
      REAL(dp)              :: rad
      INTEGER               :: mres, nres, istat, jstart
!-----------------------------------------------
!
!     xc(1,...,1) ARE THE RESONANT AMPLITUDES OF E||
!     eparmn      RESONANT RESPONSE TO xc
!     epar        IS THE RESULTING E*B/B^2 IN REAL SPACE
!
      ALLOCATE (eparmn(ns,0:mpol,-ntor:ntor), stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error in Compute_Resonant_Epar'
      eparmn = 0
      epar = 0

      DO mres=0,mpol
         DO nres=-ntor,ntor
            jstart = 2
            DO WHILE (jstart .lt. ns-1)
               CALL FindResonance(mres, nres, rad, jstart)
               IF (jstart .ge. ns-1) EXIT
               CALL AddResonantE (eparmn(:,mres,nres), rad, jstart)
            END DO
         END DO
      END DO

      DEALLOCATE(eparmn)

      END SUBROUTINE Compute_Resonant_Epar

      SUBROUTINE AddResonantE(eres, rad, jstart)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(inout) :: eres(ns)
      REAL(dp), INTENT(in)    :: rad
      INTEGER,  INTENT(in)    :: jstart
!-----------------------------------------------

      END SUBROUTINE AddResonantE


      END MODULE PERTURBATION
