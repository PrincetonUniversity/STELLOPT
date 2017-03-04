      SUBROUTINE evolve(time_step, ier_flag, liter_flag, lscreen)
      USE vmec_main
      USE vmec_params, ONLY: bad_jacobian_flag, successful_term_flag,
     1                       norm_term_flag
      USE vsvd
      USE xstuff
      USE precon2d, ONLY: ictrl_prec2d, l_comp_prec2D, compute_blocks
      USE timer_sub
      USE gmres_mod
      USE parallel_include_module
      USE vmec_params, ONLY: ntmax
!  Comment Out below JDH 2010-08-03
!      USE vmec_history
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec) :: time_step          !, r0dot
      INTEGER, INTENT(inout) :: ier_flag
      LOGICAL, INTENT(inout) :: liter_flag
      LOGICAL, INTENT(in)  :: lscreen
C-----------------------------------------------``
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      CHARACTER(LEN=*), PARAMETER :: fcn_message =
     1 "External calls to FUNCT3D: "
!      REAL(rprec), PARAMETER :: r0dot_threshold = 5.E-06_dp
      REAL(rprec) :: fsq1, dtau, b1, bprec, fac
      LOGICAL :: lfinal_mesh
      INTEGER :: lcount
      INTEGER, SAVE :: iter_on
#if defined(SKS)      
      INTEGER :: i, j, k, l, lk
      REAL(rprec) :: f3dt1, f3dt2
#endif
C-----------------------------------------------
!     IF TROUBLE CONVERGING, TRY TO RECOMPUTE PRECONDITIONER ONCE MORE...
!      IF (ictrl_prec2d.eq.1 .and. iter2.eq.(iter_on+40)) 
!     1     ictrl_prec2d = 0

!  JDH 2011-09-15 Add condition to lfinal_mesh, that iter2 - iter1 > 5
!    (5 was picked out of a hat)
!    Purpose is to keep preconditioning from being turned on immediately upon
!    a restart (V3FIT)
!   The final .and. clause is a bit complicated. The purpose of the
!   .not. lv3fit is so that if v3fita is not running, the final .and. clause
!   is always true. Read as "and, if lv3fit, then must also have iter2 - iter1 > 5"
!      lfinal_mesh = (ns .eq. ns_maxval) .and. (ictrl_prec2d.eq.0)
!     1              .and. (itype_precon.ne.0)


!      CALL second0(skston)

      lfinal_mesh = (ns .eq. ns_maxval) .and. (ictrl_prec2d.eq.0)
     1              .and. (itype_precon.ne.0)
     2              .and. ((.not. l_v3fit) .or. iter2 - iter1 .ge. 5)

      IF (iter2 .lt. 10) THEN
         ictrl_prec2d = 0
         lqmr = .false.
         iter_on = -1
      ELSE IF (lfinal_mesh .and. 
     1        (fsqr+fsqz+fsql).lt.prec2d_threshold) THEN
!     2                     .and. r0dot.lt.r0dot_threshold) THEN
         lqmr = (itype_precon .ge. 2)
         lfirst = (lqmr .and. iter_on.eq.-1) 

!
!        INITIATES 2D PRECONDITIONER CALCULATION
!
         IF (iter_on .eq. -1) THEN
            IF (lqmr) THEN
               nstep = 5
               niter = iter2+100                   !Limit to 100 preconditioner steps
            END IF
            iter_on = iter2                     !Flag to monitor progress of preconditioner
         ELSE
            iter_on = iter2-11
         END IF

!         iter_on = iter2                     !Flag to monitor progress of preconditioner

!SPH022111: ADD NEW CONTROL PARAMETER, l_comp_prec2D, TO FORCE RECALCULATION
!           OF PRECONDITIONING BLOCKS IN V3FIT, FOR EXAMPLE
         IF (lfirst .or. l_comp_prec2D) THEN
            IF (l_v3fit) WRITE(*,*) 'VMEC Evolve:compute_blocks'
            CALL compute_blocks (xc,xcdot,gc)
         ENDIF
         IF(l_v3fit) WRITE(*,*) 'VMEC Evolve:prec2d_On iter2 =', iter2
         l_comp_prec2D = .FALSE.
         ictrl_prec2d = 1
         time_step = 0.50_dp
         iter1 = iter2-1; fsq = fsqr1 + fsqz1 + fsql1
         IF(PARVMEC) THEN
           pxstore = pxc
           pxcdot = 0
         ELSE
           xstore = xc
           xcdot = 0
         END IF
      END IF

!
!     COMPUTE MHD FORCES
!     MUST CALL funct3d EVEN WHEN IN 2D PRECONDITIONING MODE, SINCE
!     INITIAL RESIDUALS MUST BE KNOWN WHEN CALLING gmres_fun, etc.
!
#if defined (SKS)      
      CALL second0(f3dt1)
      f3d_num(NS_RESLTN) = f3d_num(NS_RESLTN)+1
#endif
      IF (PARVMEC) THEN
        CALL funct3d_par (lscreen, ier_flag)
      ELSE
        CALL funct3d (lscreen, ier_flag)
      END IF
#if defined (SKS)      
      CALL second0(f3dt2)
      f3d_time(NS_RESLTN) = f3d_time(NS_RESLTN) + (f3dt2-f3dt1)
      funct3d_time = funct3d_time + (f3dt2-f3dt1)
#endif        

!
!     COMPUTE ABSOLUTE STOPPING CRITERION
      IF (iter2.eq.1 .and. irst.eq.2) THEN
         ier_flag = bad_jacobian_flag
!  JDH 2012-04-24. Revise this absolute stopping criterion, so that if v3fit
!    is running, then have to iterate at least 2 * nvacskip steps
!    (2 picked out of a hat) (nvacskip - to make sure vacuum gets updated)
!    before returning.
!      ELSE IF (fsqr.le.ftolv .and. fsqz.le.ftolv .and.
!     1         fsql.le.ftolv) THEN
      ELSE IF (fsqr.le.ftolv .and. fsqz.le.ftolv .and.
     1         fsql.le.ftolv .and.
     2         ((.not. l_v3fit) .or. iter2 - iter1 .ge. 2 * nvacskip))
     3         THEN
         liter_flag = .false.
         ier_flag = successful_term_flag
         IF (lqmr) THEN
            WRITE (nthreed,'(/,2x,a,i5)') fcn_message,nfcn
            WRITE (*,'(/,2x,a,i5)') fcn_message,nfcn
         END IF
      ENDIF

      IF (ier_flag.ne.norm_term_flag .or. .not.liter_flag) RETURN

      IF (lqmr) THEN
         IF (irst .eq. 2) THEN
!           TRAP IRST=2, CHANGE IN GSQRT SIGN, IF INITIAL PERT TOO LARGE
            xc = xsave + 0.25_dp*(xc - xsave)
            RETURN
         END IF

         CALL gmres_fun(ier_flag, itype_precon-1)

!  JDH 2011-06-20. Comment out V3FIT stabilization
!V3FITA STABILIZATION (021711): NOTE - lower 100 for turning off precondition sooner
C         b1 = (fsqr+fsqz+fsql)/prec2d_threshold
C         IF (b1 .ge. 10) THEN
C            xc = xstore + (xc-xstore)/b1
C            IF (b1 .gt. 1000) THEN
C               xc = xstore
C               lqmr = .FALSE.
C               ictrl_prec2d = 0
C               time_step = 0.9*time_step
C               RETURN
C            END IF
C         END IF
!END OF V3FITA STABILIZATION

         IF (lfreeb) RETURN
         DO lcount = ns, 2*irzloff, ns
            IF (xsave(lcount) .ne. xc(lcount))
     1         PRINT *,' xsave = ',xsave(lcount),' != xc = ',
     2         xc(lcount),' for lcount = ',lcount
         END DO
         RETURN
      END IF

!     COMPUTE DAMPING PARAMETER (DTAU) AND EVOLVE
!     R, Z, AND LAMBDA ARRAYS IN FOURIER SPACE

      fsq1 = fsqr1 + fsqz1 + fsql1

      IF (iter2 .eq. iter1) otau(:ndamp) = cp15/time_step

      IF (ictrl_prec2d .eq. 0) THEN
         bprec = 1
      ELSE
         bprec = 3
      END IF

      IF (iter2.gt.iter1 .and. fsq1.ne.zero) 
     1    dtau = MIN(ABS(LOG(fsq1/fsq)), bprec*cp15)

      fsq = fsq1

      otau(1:ndamp-1) = otau(2:ndamp)

      IF (iter2 .gt. iter1) otau(ndamp) = dtau/time_step
!REMOVED 071505: OTHERWISE I=1 STATE REPEATED (SKIP THIS TO GET OUT OF ITER2=1 STATE)
!     IF (iter2 .le. 1) RETURN

      otav = SUM(otau(:ndamp))/ndamp
      dtau = time_step*otav/2

      b1  = one - dtau
      fac = one/(one + dtau)

!
!     THIS IS THE TIME-STEP ALGORITHM. IT IS ESSENTIALLY A CONJUGATE
!     GRADIENT METHOD, WITHOUT THE LINE SEARCHES (FLETCHER-REEVES),
!     BASED ON A METHOD GIVEN BY P. GARABEDIAN

!

#if defined (SKS)
      IF(PARVMEC) THEN
        IF (lactive) THEN 
          lk=0
          DO l=1, 3*ntmax
            DO i=t2lglob, t2rglob
              DO k=0, mpol1
                DO j=0, ntor
                  lk = j + (ntor+1)*k +
     1           (ntor+1)*(mpol1+1)*(i-1)+(ntor+1)*(mpol1+1)*ns*(l-1)+1
                  pxcdot(lk) = fac*(b1*pxcdot(lk) + time_step*pgc(lk))
                  pxc(lk)    = pxc(lk) + time_step*pxcdot(lk)
                END DO
              END DO
            END DO
          END DO
        END IF
      ELSE
#endif
        xcdot = fac*(b1*xcdot + time_step*gc)
        xc    = xc + time_step*xcdot
#if defined (SKS)
      ENDIF
#endif
      END SUBROUTINE evolve
