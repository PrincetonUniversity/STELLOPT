
SUBROUTINE flint_bo()
! eps_eff Integration
! and additional output specified in the paper
!
! Input variables:  bmref    -  $B_0$ - reference value of mod-B
!                   rt0      -  $R$   - reference value of big radius
!                   phi0     -  starting value of $\phi$
!                   theta0   -  starting value of $\theta$
!
! Output variables: epspar(i)-  partial contributions to effective ripple
!                               from different classes of trapped particles
!                               index i=1,...,multra  means single-trapped,
!                               double-trapped,...,up to multra-1,
!                               epspar(multra) contains the contributions
!                               from all upper classes
!                   epstot   -  total effective ripple - $\epsilon_{eff}^{3/2}$,
!                               the main result of computation
!                   ctrone   -  fraction of particles trapped once (see Eq.(36))
!                   ctrtot   -  fraction of all trapped particles  (see Eq.(36))
!                   bareph   -  $\bar \epsilon_h$ - 'ripple' amplitude defined
!                               through fraction of single-trapped particles
!                   barept   -  $\bar \epsilon_t$ - 'toroidal' amplitude
!                               defined through full fraction of trapped prts
!
! **********************************************************************
! Modules
! **********************************************************************
  USE neo_precision
  USE neo_units
  USE neo_exchange
  USE neo_output
  USE partpa_bo
  USE neo_rk4dbo
  USE safe_open_mod                                   ! SPH
! **********************************************************************
! Local definitions
! **********************************************************************
  IMPLICIT NONE
!
  REAL(kind=dp) ::   etamin,etamax,hphi
  REAL(kind=dp) ::   phi
  REAL(kind=dp) ::   heta
  REAL(kind=dp) ::   coeps
  REAL(kind=dp) ::   adimax,aditot,adimax_s,aditot_s
  REAL(kind=dp) ::   y2,y3
  REAL(kind=dp) ::   y2_s,y3_s,y4,y4_s,y3npart,y3npart_s
  REAL(kind=dp) ::   bmod,gval,geodcu,qval
  REAL(kind=dp) ::   acc_d, acc_m
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE, TARGET ::  y
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE ::  bigint, bigint_s
  INTEGER                                  ::  i,m_cl,j1,n, istat
  INTEGER,       DIMENSION(:), ALLOCATABLE ::  iswst

  REAL(kind=dp)                            ::  epstot_check
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE ::  epspar_check

! INTEGER,       PARAMETER              ::  ntheta_bins = 100
! INTEGER                               ::  theta_bind
! INTEGER,       DIMENSION(ntheta_bins) ::  theta_bins
! INTEGER,       DIMENSION(ntheta_bins) ::  theta_bins_c
! REAL(kind=dp)                         ::  delta_theta
  REAL(kind=dp)                         ::  theta
!
  REAL(kind=dp)                         ::  theta0, phi0
  REAL(kind=dp)                         ::  theta_rs, theta_iota
  REAL(kind=dp)                         ::  theta_d,theta_d_min
  REAL(kind=dp)                         ::  theta_gap
  REAL(kind=dp)                         ::  iota_bar_fp

  INTEGER                               ::  istepc
  INTEGER                               ::  m, m_iota, n_iota
  INTEGER                               ::  nfl
  INTEGER                               ::  n_gap
  INTEGER                               ::  nstep_max_c
  INTEGER                               ::  exist_first_ratfl
  INTEGER                               ::  max_class
  REAL(kind=dp)                         ::  j_iota_i
  REAL(kind=dp)                         ::  add_on

  REAL(kind=dp),               POINTER  ::  p_theta, p_bm2
  REAL(kind=dp),               POINTER  ::  p_bm2gv, p_bm3ge
  REAL(kind=dp), DIMENSION(:), POINTER  ::  p_i, p_h
! **********************************************************************
! Allocation
! **********************************************************************
  ndim = npq + 2 * npart
  ALLOCATE( isw(npart) )
  ALLOCATE( ipa(npart) )
  ALLOCATE( icount(npart) )
  ALLOCATE( iswst(npart) )
  ALLOCATE( eta(npart) )

  ALLOCATE( y(ndim) )

  IF (.NOT. ALLOCATED(epspar) ) ALLOCATE( epspar(multra) )
  ALLOCATE( bigint(multra) )
  ALLOCATE( bigint_s(multra) )
  ALLOCATE( epspar_check(multra) )
! **********************************************************************
! Pointer
! **********************************************************************
  p_theta => y(1)
  p_bm2   => y(2)
  p_bm2gv => y(3)
  p_bm3ge => y(4)
  p_i     => y(npq+1:npq+npart)
  p_h     => y(npq+npart+1:npq+2*npart)
! **********************************************************************
! Initial settings
! **********************************************************************
  ierr   = 0
  hphi   = twopi/(nstep_per*nper)
  bmod0  = bmref
!
  exist_first_ratfl = 0
  hit_rat = 0
!
  theta0 = theta_bmax
  phi0   = phi_bmax
! **********************************************************************
! Definition of eta-values for particles

! **********************************************************************
  etamin = b_min / bmref
  etamax = b_max / bmref
  heta=(etamax-etamin)/(npart-1)
!  heta=(etamax-etamin)/(npart+1)
!  etamin = etamin + heta
!  etamax = etamax - heta
  etamin = etamin + heta / 2
  DO i=1,npart
     eta(i)=etamin+heta*(i-1)
  ENDDO
  coeps=pi*rt0**2*heta/(8*sqrt(2._dp))
  j_iota_i = curr_pol(psi_ind) + iota(psi_ind) * curr_tor(psi_ind)
! **********************************************************************
  isw    = 0
  ipa    = 0
  iswst  = 0
  ipmax  = 0
  icount = 0
  bigint = 0
  adimax = 0
  aditot = 0
! **********************************************************************
  nfp_rat = 0
  nfl_rat = 0
! **********************************************************************
  phi     = phi0
  y       = 0
  p_theta = theta0
! **********************************************************************
! Bins for theta at periodicity boundary
! **********************************************************************
! delta_theta = twopi / no_bins
! theta_bins  = 0
! theta_bind  = MOD(FLOOR(theta0/delta_theta),no_bins)+1
! theta_bins(theta_bind) = 1
!
  theta_d_min = twopi
! **********************************************************************
  IF (write_integrate == 1) THEN
     OPEN(unit=w_u5,file='conver.dat',status='replace',form='formatted')
  ENDIF
  IF (write_diagnostic == 1) THEN
     OPEN(unit=w_u7,file='diagnostic.dat',status='replace',form='formatted')
  ENDIF
! **********************************************************************
! Evaluation of Splines
! **********************************************************************
  CALL neo_eval(theta0,phi0,bmod,gval,geodcu,pard0,qval)
! **********************************************************************
! Intergration steps (Summation)
! **********************************************************************
  nstep_max_c = nstep_max
  istepc = 0
  max_class = 0
  DO n=1,nstep_max
     DO j1=1,nstep_per
        CALL RK4D_BO(phi,y,hphi)
        DO i=1,npart
           IF(isw(i).EQ.2) THEN
              m_cl = MIN(multra,ipa(i))
              m_cl = MAX(1,m_cl)
              IF(ipa(i).EQ.1) adimax = p_i(i)
              add_on = p_h(i)**2 / p_i(i) * iswst(i)
              bigint(m_cl) = bigint(m_cl) + add_on

              IF ( write_diagnostic .EQ. 1 .AND. iswst(i) .EQ. 1) THEN
                 istepc = istepc + 1
                 WRITE (w_u7,'(1x,i8,1x,i8,1x,i8,1x,e20.10)')          &
                      i,icount(i),ipa(i),add_on
                 IF ( ipa(i) .GT. max_class ) max_class = ipa(i)
              END IF

              iswst(i) = 1
              p_h(i) = 0
              p_i(i) = 0
              isw(i) = 0
              icount(i) = 0
              ipa(i) = 0
           END IF
        END DO
!
        IF(ipmax.EQ.1) THEN
           ipmax=0
           aditot=aditot+adimax
        END IF
     END DO
! **********************************************************************
! Checks after each field-period
! **********************************************************************
     IF (write_integrate == 1) THEN
        epstot_check = 0
        DO m_cl=1,multra
           epspar_check(m_cl) = coeps*bigint(m_cl)*p_bm2/p_bm2gv**2
           epstot_check = epstot_check + epspar_check(m_cl)
        ENDDO
        WRITE(w_u5,format220)                                          &
             DBLE(n), epstot_check, p_bm3ge,                           &
             p_i(npart)/p_bm2, aditot/p_bm2
     END IF
! **********************************************************************
! Bins for theta at periodicity boundary
! **********************************************************************
     theta = p_theta
!    theta_bind = MOD(FLOOR(theta/delta_theta),ntheta_bins)+1
!    theta_bins(theta_bind) = theta_bins(theta_bind) + 1
!    acc_d = (MAXVAL(theta_bins)-MINVAL(theta_bins)) *                 &
!         SQRT(DBLE(ntheta_bins)) / n / 2
!    theta_bins_c = 0
!    WHERE (theta_bins .NE. 0)
!       theta_bins_c = 1
!    END WHERE
!    acc_m = (DBLE(SUM(theta_bins))/DBLE(SUM(theta_bins_c)) -          &
!             MINVAL(theta_bins)) *                                    &
!             SQRT(DBLE(ntheta_bins)) / n
! **********************************************************************
! Check for (near) rational surfaces
! **********************************************************************
     IF (n .LE. nstep_min) THEN
        theta_rs = theta - theta0
        IF (n .EQ. 1) THEN
           theta_iota = theta_rs
           iota_bar_fp = theta_iota / twopi
        END IF
        m = FLOOR(theta_rs / twopi)
        theta_rs = theta_rs - m * twopi
        IF (theta_rs .LE. pi) THEN
           theta_d = theta_rs
        ELSE
           theta_d = theta_rs - twopi
        END IF
        IF (ABS(theta_d) .LT. ABS(theta_d_min)) THEN
           theta_d_min = theta_d
           n_iota = n
           IF (theta_d .GE. ZERO) THEN
              m_iota = m
           ELSE
              m_iota = m+1
           END IF
        END IF
     END IF
     IF (n .EQ. nstep_min) THEN
        theta_gap = twopi / n_iota
        n_gap = n_iota * INT(ABS(theta_gap / theta_d_min))
        IF (n_gap .GT. nstep_min) THEN
           nstep_max_c = n_gap
        ELSE
           nstep_max_c = n_gap * CEILING( real(nstep_min,kind=dp) / n_gap )
        END IF

!!$        PRINT *, 'theta_iota:      ',theta_iota
!!$        PRINT *, 'theta_rs:        ',theta_rs
!!$        PRINT *, 'theta_d:         ',theta_d
!!$        PRINT *, 'theta_d_min:     ',theta_d_min
!!$        PRINT *, 'theta_gap:       ',theta_gap
!!$        PRINT *, 'n_iota:          ',n_iota
!!$        PRINT *, 'n_gap:           ',n_gap
!!$        PRINT *, 'nstep_min:       ',nstep_min
!!$        PRINT *, 'n_gap:           ',n_gap
!!$        PRINT *, 'DBLE(nstep_min): ',DBLE(nstep_min)
!!$        PRINT *, 'DBLE(n_gap):     ',DBLE(n_gap)
!!$        PRINT *, 'DBLE/DBLE:       ',DBLE(nstep_min)/DBLE(n_gap)
!!$        PRINT *, 'CEILING:         ',CEILING(DBLE(nstep_min) / DBLE(n_gap))
!!$        PRINT *, 'nstep_max_c:     ',nstep_max_c

! **********************************************************************
! Decision about Rational Surfaces
! **********************************************************************
        IF (nstep_max_c .GT. nstep_max) THEN
           hit_rat = 1
           nfp_rat = CEILING(ONE / acc_req / iota_bar_fp)
           IF (MODULO(nfp_rat,n_iota) .NE. 0) THEN
              nfp_rat = nfp_rat + n_iota - MODULO(nfp_rat,n_iota)
           END IF
           IF (nfp_rat .GE. nstep_min) THEN
              exist_first_ratfl = 1
              nstep_max_c = nfp_rat
           END IF
           nfl_rat = CEILING(real(no_bins,kind=dp) / n_iota)
           delta_theta_rat =  theta_gap / (nfl_rat + 1)
           IF ( calc_nstep_max .EQ. 1 ) hit_rat = 0
           IF ( hit_rat .EQ. 1 .AND. exist_first_ratfl .EQ. 0 ) EXIT
        END IF
     END IF
! **********************************************************************
! Accuracy
! **********************************************************************
!     IF (acc_d .LT. acc_req .AND. n .GT. nstep_min) THEN
!        EXIT
!     END IF
     IF ( calc_nstep_max .EQ. 0 .AND. n .EQ. nstep_max_c ) EXIT
  END DO
  nintfp = n
!
  y2      = p_bm2
  y3      = p_bm2gv
  y4      = p_bm3ge
  y3npart = p_i(npart)
! **********************************************************************
  IF (write_integrate == 1) THEN
     CLOSE(unit=w_u5)
  END IF
! **********************************************************************
! Calculations for rational surfaces
! **********************************************************************
  IF (hit_rat .EQ. 1) THEN

     IF (exist_first_ratfl .EQ. 0) THEN
        IF (write_integrate == 1) THEN
           OPEN(unit=w_u5,file='conver.dat',status='replace',form='formatted')
        ENDIF
        bigint  = 0
        adimax  = 0
        aditot  = 0
        y2      = 0
        y3      = 0
        y4      = 0
        y3npart = 0
     ELSE
        IF (write_integrate == 1) THEN
           OPEN(unit=w_u5,file='conver.dat',status='old',                   &
                position='append',form='formatted')
        END IF
     END IF
     DO nfl = exist_first_ratfl,nfl_rat
        bigint_s  = 0
        adimax_s  = 0
        aditot_s  = 0
        isw       = 0
        ipa       = 0
        iswst     = 0
        ipmax     = 0
        icount    = 0
        phi       = phi0
        y         = 0
        theta     = theta0 + nfl * delta_theta_rat
        p_theta   = theta
        DO n=1,nfp_rat
           CALL neo_eval(theta,phi,bmod,gval,geodcu,pard0,qval)
           DO j1=1,nstep_per
              CALL RK4D_BO(phi,y,hphi)
              DO i=1,npart
                 IF(isw(i).EQ.2) THEN
                    isw(i) = 0
                    icount(i) = 0
                    m_cl = MIN(multra,ipa(i))
                    m_cl = MAX(1,m_cl)
                    IF(ipa(i).EQ.1) adimax_s = p_i(i)
                    ipa(i) = 0
                    bigint_s(m_cl) = bigint_s(m_cl)+p_h(i)**2/p_i(i)*iswst(i)
                    iswst(i) = 1
                    p_h(i) = 0
                    p_i(i) = 0
                 END IF
              END DO
              IF(ipmax.EQ.1) THEN
                 ipmax    = 0
                 aditot_s = aditot_s + adimax_s
              END IF
           END DO

           IF (write_integrate == 1) THEN
              epstot_check = 0
              DO m_cl=1,multra
                 epspar_check(m_cl) = coeps*bigint_s(m_cl)*p_bm2/p_bm2gv**2
                 epstot_check = epstot_check + epspar_check(m_cl)
              ENDDO
              WRITE(w_u5,format220)                                         &
                   DBLE(nfl*nfp_rat + n),                                   &
                   epstot_check,                                            &
                   p_bm3ge,p_i(npart)/p_bm2,aditot_s/p_bm2
           END IF

        END DO

        y2_s      = p_bm2
        y3_s      = p_bm2gv
        y4_s      = p_bm3ge
        y3npart_s = p_i(npart)

        bigint  = bigint + bigint_s
        aditot  = aditot + aditot_s
        y2      = y2 + y2_s
        y3      = y3 + y3_s
        y4      = y4 + y4_s
        y3npart = y3npart + y3npart_s
     END DO
     n = nfp_rat * (nfl_rat + 1)
  END IF
  IF (write_integrate == 1) THEN
     CLOSE(unit=w_u5)
  END IF
  IF (write_diagnostic == 1) THEN
     CLOSE(unit=w_u7)
     OPEN(unit=w_u7,file='diagnostic_add.dat',status='replace',form='formatted')
     WRITE(w_u7,'(4(1x,i8),6(1x,e20.10))')      &
          psi_ind,istepc,npart,max_class,b_min,b_max,bmref,coeps,y2,y3
     CLOSE(unit=w_u7)
  END IF
! **********************************************************************
! Final results
! **********************************************************************
  epstot=0
  DO m_cl=1,multra
     epspar(m_cl) = coeps*bigint(m_cl)*y2/y3**2
     epstot = epstot + epspar(m_cl)
  ENDDO
!
  ctrone = aditot / y2
  ctrtot = y3npart / y2
!
  bareph=(pi*ctrone)**2/8
  barept=(pi*ctrtot)**2/8
!
  drdpsi=y2/y3
!
  yps = y4 * j_iota_i
! **********************************************************************
! Log File
! **********************************************************************
  IF (w_u6_open .EQ. 0) THEN
     call safe_open(w_u6, istat, 'neolog.' // trim(extension), 'replace', 'formatted')
     w_u6_open = 1
  ELSE
     OPEN(unit=w_u6,file='neolog.' // trim(extension), status='old',   &
          iostat=istat,position='append',form='formatted')
     if (istat .ne. 0) stop 'NEOLOG.DAT CANNOT BE OPENED IN NEO FLINT_BO'
  END IF
  WRITE(w_u6,'(7(1x,i6),1x,1(d16.8))')                                 &
       psi_ind,n_iota,m_iota,n_gap,nfp_rat,nfl_rat+1,n,epstot
  CLOSE(w_u6)
! **********************************************************************
! Deallocation
! **********************************************************************
  DEALLOCATE( isw )
  DEALLOCATE( ipa )
  DEALLOCATE( icount )
  DEALLOCATE( iswst )
  DEALLOCATE( eta )

  DEALLOCATE( y )

  DEALLOCATE( bigint )
  DEALLOCATE( bigint_s )
  DEALLOCATE( epspar_check )

  RETURN
END SUBROUTINE flint_bo
