      SUBROUTINE chisq_bmn_tgt (target, sigma, n_tgt, m_tgt, n_n,
     1                           hs, num, mnboz, nsval, nopt, lscreen)
      USE stel_kinds
      USE chisq_mod
      USE optim, ONLY: nfp_opt, bigno
      USE optim_params, ONLY: helicity
      USE read_boozer_mod, ONLY: bmn_b=>bmnc_b, ixm_b, ixn_b
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER,     INTENT(in)    :: mnboz, nsval, nopt,
     1                              n_n, n_tgt(n_n), m_tgt(n_n)
      INTEGER,     INTENT(inout) :: num
      REAL(rprec), INTENT(in)    :: hs, sigma(n_n), target(n_n)
      LOGICAL,     INTENT(in)    :: lscreen 
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, one = 1, c1p5 = 1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: jcount, m, n, mn, num0, i, j
      REAL(rprec) :: bmax, sj
C-----------------------------------------------
!
!     Computes the ratio of energy (B**2) in modes with the undesirable helicities (match)
!     to the energy in modes with the desired helicity (bnorm) and minimizes that...
!     PORTED FROM PPPL Code -SAL 08/04/11
      

      IF (nopt .gt. 0) THEN
         
         bmax  = MAXVAL(ABS(bmn_b(1:mnboz,nsval)))
         
         T_LOOP: DO i = 1, n_n
         
            IF ( sigma(i) .ge. bigno ) CYCLE T_LOOP
            
            DO mn = 1, mnboz
               
               IF ( ixn_b(mn)/nfp_opt == n_tgt(i) .and.
     1              ixm_b(mn)         == m_tgt(i) ) THEN
         
                  num               = num + 1
                  wegt(num)         = sigma(i)
                  index_array       = ivar_bmn
                  chisq_target(num) = target(i)
                  chisq_match(num)  = bmn_b(mn,nsval) / bmax
                  
                  CYCLE T_LOOP
                  
               END IF
               
            END DO
            
            ! If there are no matches
            IF ( lscreen ) THEN
               
               WRITE(*,*)
     1           ' Chisq_bmn_tgt: harmonic targeted but not calculated!'
               WRITE(*,*) ' n = ', n_tgt(i), '  m = ', m_tgt(i)
               
            END IF
            
            num               = num + 1
            wegt(num)         = sigma(i)
            index_array(num)  = ivar_bmn
            chisq_target(num) = target(i)
            chisq_match(num)  = 0
            
         END DO T_LOOP

      ELSE
         num0 = COUNT( sigma(:) < bigno)
         IF (nopt .eq. -2) 
     1             chisq_descript(num+1:num+num0) = descript(ivar_bmn)
         num = num + num0
      END IF


      END SUBROUTINE chisq_bmn_tgt
