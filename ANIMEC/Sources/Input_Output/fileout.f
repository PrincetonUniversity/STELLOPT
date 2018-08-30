      SUBROUTINE fileout(iseq, ictrl_flag, ier_flag, lscreen)
      USE vmec_main
      USE vac_persistent
      USE realspace
      USE vmec_params, ONLY: mscale, nscale, signgs, uminus,
     1      norm_term_flag, more_iter_flag, output_flag, 
     2      cleanup_flag, successful_term_flag
      USE vforces
      USE vsvd
      USE xstuff, ONLY: xc, gc, xsave, scalxc
      USE precon2d, ONLY: ictrl_prec2d
      USE timer_sub
#ifdef _HBANGLE
      USE angle_constraints, ONLY: free_multipliers, getrz
#endif
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: iseq, ictrl_flag
      INTEGER, INTENT(inout) :: ier_flag
      LOGICAL :: lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER :: istat, loc_ier_flag
      LOGICAL, PARAMETER :: lreset_xc = .false.
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: js, istat1=0, irst0
      REAL(rprec), DIMENSION(:), POINTER :: lu, lv
      REAL(rprec), ALLOCATABLE :: br_out(:), bz_out(:)
      CHARACTER(LEN=*), PARAMETER, DIMENSION(0:10) :: werror = (/
     1   'EXECUTION TERMINATED NORMALLY                            ',
     2   'INITIAL JACOBIAN CHANGED SIGN (IMPROVE INITIAL GUESS)    ',
     3   'FORCE RESIDUALS EXCEED FTOL: MORE ITERATIONS REQUIRED    ',
     4   'VMEC INDATA ERROR: NCURR.ne.1 but BLOAT.ne.1.            ',
     5   'MORE THAN 75 JACOBIAN ITERATIONS (DECREASE DELT)         ',
     6   'ERROR READING INPUT FILE OR NAMELIST                     ',
     7   'NEW AXIS GUESS STILL FAILED TO GIVE GOOD JACOBIAN        ',
     8   'PHIEDGE HAS WRONG SIGN IN VACUUM SUBROUTINE              ',
     9   'NS ARRAY MUST NOT BE ALL ZEROES                          ',
     A   'ERROR READING MGRID FILE                                 ',
     B   'VAC-VMEC I_TOR MISMATCH : BOUNDARY MAY ENCLOSE EXT. COIL ' /)
      CHARACTER(LEN=*), PARAMETER ::
     1    Warning = " Error deallocating global memory FILEOUT"
      LOGICAL :: log_open, lwrite, loutput, lterm
C-----------------------------------------------
      lu => czmn;   lv => crmn

!
!     COMPUTE REMAINING COVARIANT COMPONENT OF B (BSUBS),
!     CYLINDRICAL COMPONENTS OF B (BR, BPHI, BZ), AND
!     AVERAGE EQUILIBRIUM PROPERTIES AT END OF RUN
!
      iequi = 1
      lterm = ier_flag.eq.norm_term_flag  .or.
     1        ier_flag.eq.successful_term_flag
      lwrite = lterm .or. ier_flag.eq.more_iter_flag
      loutput = (IAND(ictrl_flag, output_flag) .ne. 0)
      loc_ier_flag = ier_flag
      if (ier_flag .eq. successful_term_flag) 
     1    loc_ier_flag = norm_term_flag

      IF (lwrite .and. loutput) THEN
!
!     The sign of the jacobian MUST multiply phi to get the physically
!     correct toroidal flux
!
         phi(1) = zero
         DO js = 2, ns
            phi(js) = phi(js-1) + phip(js)
         END DO
         phi = (signgs*twopi*hs)*phi

!        Must save irst value if in "restart" mode
         irst0 = irst
         CALL funct3d (lscreen, istat)

!        Write out any special files here
!        CALL dump_special

         irst = irst0

         CALL second0 (teqfon)
         ALLOCATE(br_out(nrzt), bz_out(nrzt), stat=istat)
         gc = xc
#ifdef _HBANGLE
         CALL getrz(gc)
#endif
         CALL eqfor (br_out, bz_out, clmn, blmn, rcon(1,1), 
     1               gc, ier_flag)
         CALL second0 (teqfoff)
         timer(teqf) = timer(teqf) + teqfoff - teqfon
      END IF
  
!      CALL free_mem_precon
!
!     Call WROUT to write output or error message if lwrite = false
!

      IF (loutput .and. ASSOCIATED(bzmn_o)) THEN
         CALL second0 (twouton)

         CALL wrout (bzmn_o, azmn_o, clmn, blmn, crmn_o, czmn_e,
     1        crmn_e, xsave, gc, loc_ier_flag, lwrite
#ifdef _ANIMEC
     2       ,brmn_o, sigma_an, ppar, pperp, onembc, pp1, pp2, pp3
#endif 
     3        )
         CALL second0 (twoutoff)

         timer(twout) = timer(twout) + twoutoff - twouton

         IF (ntor .eq. 0) THEN
            CALL write_dcon (xc)
         END IF

         IF (lscreen .and. ier_flag.ne.more_iter_flag) 
     1   PRINT 120, TRIM(werror(loc_ier_flag))
         IF (lscreen .and. lterm) 
     1   PRINT 10, TRIM(input_extension), ijacob

         IF (nthreed .gt. 0) THEN
            WRITE (nthreed,120) TRIM(werror(loc_ier_flag))
            IF (.not. lterm) GOTO 1000
            WRITE (nthreed, 10) TRIM(input_extension), ijacob
            CALL write_times(nthreed, lscreen, lfreeb, lrecon, 
     1                       ictrl_prec2d.ne.0)
         END IF
      END IF

   10 FORMAT(' FILE : ',a,/,' NUMBER OF JACOBIAN RESETS = ',i4,/)
  120 FORMAT(/1x,a,/)

!
!     TESTING READ_WOUT MODULE WRITING ROUTINES
!
      IF (ALLOCATED(br_out)) THEN
!        IF (lscreen) CALL TestWout(xc, br_out, bz_out, crmn_e, czmn_e)
         DEALLOCATE (br_out, bz_out)
      END IF

!     END TEST

!
!     WRITE SEQUENCE HISTORY FILE
!
      INQUIRE(unit=nlog,opened=log_open)
      IF (lrecon .and. log_open) THEN
         IF (iseq .eq. 0) WRITE (nlog, 100)
         WRITE (nlog, 110) iseq + 1, iter2, total_chi_square_n,
     1      1.e-6_dp*ctor/mu0, 1.e-3_dp*ppeak/mu0, torflux, r00,
     2      timer(tsum), input_extension
      ENDIF

  100 FORMAT(' SEQ ITERS  CHISQ/N',
     1   '  TORCUR  PRESMAX  PHIEDGE     R00 CPU-TIME  EXTENSION')
  110 FORMAT(i4,i6,f8.2,3f9.2,f8.2,f9.2,2x,a20)

 1000 CONTINUE

!
!     DEALLOCATE GLOBAL MEMORY AND CLOSE FILES
!
      IF (IAND(ictrl_flag, cleanup_flag).eq.0 .or. 
     1         ier_flag.eq.more_iter_flag) RETURN

      IF (ALLOCATED(cosmu))
     1  DEALLOCATE(cosmu, sinmu, cosmum, sinmum, cosmui, cosmumi,
     2  sinmui, sinmumi, cosnv, sinnv, cosnvn, sinnvn, cosmui3,
     3  cosmumi3, cos01, sin01, stat=istat1)
      IF (istat1 .ne. 0) PRINT *, Warning // "#1"
#ifdef _HBANGLE
        CALL free_multipliers
#endif

      IF (ALLOCATED(xm)) DEALLOCATE (xm, xn, ixm, xm_nyq, xn_nyq, 
     1   jmin3, mscale, nscale, uminus, stat=istat1)
      IF (istat1 .ne. 0) PRINT *, Warning // "#2"

      IF (ALLOCATED(tanu))
     1  DEALLOCATE(tanu, tanv, sinper, cosper, sinuv, cosuv,
     2  sinu, cosu, sinv, cosv, sinui, cosui, cmns, csign, sinu1,
     3  cosu1, sinv1, cosv1, imirr, xmpot, xnpot, stat=istat1)
      IF (istat1 .ne. 0) PRINT *, Warning // "#3"

      CALL free_mem_funct3d
      CALL free_mem_ns (lreset_xc)
      CALL free_mem_nunv
      CALL free_persistent_mem

      CALL close_all_files

      END SUBROUTINE fileout
