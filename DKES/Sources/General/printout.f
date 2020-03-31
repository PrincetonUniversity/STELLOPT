      SUBROUTINE dkes_printout(fz1p, fz1m, fz3p, fz3m, srces)

c  This subroutine calculates the diffusion coefficients, parallel
c  viscous stress, and banana-plateau flux, and writes them to dkesout.
c
c      fz1p:      distribution function in response to source type = 1,+ (sin for stell sym)
c      fz1m:      distribution function in response to source type = 1,- (cos for stell sym)
c      fz3p:      distribution function in response to source type = 3,+ (sin for stell sym)
c      fz3m:      distribution function in response to source type = 3,- (cos for stell sym)
c
c      Uses Eqs (24) and (33) in W. I. van Rij, et. al. paper
c
c      See subroutine RESIDUE for definitions of residuals
c
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE Vimatrix
      USE Vnamecl2
      USE dkes_input, ONLY: ipmb, lalpha, lscreen, itab_out,     !record file addition
     1  psip, chip, btheta, bzeta                                !record file addition
      USE dkes_realspace, ONLY: bstrs, mpnt, mpnt2, mpnt3, mpnt4,
     1    diagl, DKES_rad_dex, DKES_L11p, DKES_L33p, DKES_L31p,
     2    DKES_L11m, DKES_L33m, DKES_L31m,
     3    DKES_scal11, DKES_scal33, DKES_scal31
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(mpnt,0:lalpha), INTENT(in) :: 
     &             fz1p, fz1m, fz3p, fz3m
      REAL(rprec), DIMENSION(mpnt,4,2,2), INTENT(in)    :: srces
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: l0 = 1, l1 = 2, l2 = 3, l3 = 4,
     1   pgrad = 1, epar = 2, plus = 1, minus = 2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: mn
      REAL(rprec) :: L11p, L33p, rs11p, rs33p, L13p, L31p,
     1   rs13p, rs31p,
     1   stress1, bpflux1, L11m, L33m, rs11m, rs33m,
     2   L13m, L31m, rs13m, rs31m, g31max, g31min,
     3   del, av, scal11, scal33, scal13, srcp, srcm, t1, t2, t3, t4
      REAL(kind=rprec), DIMENSION(16) :: rsds         !record file addition
      REAL(kind=rprec) :: rsds_max                    !record file addition
      INTEGER :: ir                                   !record file addition
      CHARACTER :: tb*1                               !record file addition
C-----------------------------------------------------------------------
      tb = char(9)                                   !record file addition

      srcp = SUM(fz1p(1:mpnt,l0:l0+3)*srces(1:mpnt,l0:l0+3,1,plus))
      t1   =-SUM(fz1m(1:mpnt,l0:l0+3)*
     1           srces(1:mpnt,l0:l0+3,pgrad,minus))/vp
      t2   = SUM(fz1p(1:mpnt,l0)*srces(1:mpnt,l0,1,plus))/vp
      t3   =-s1cs1
      t4   = SUM(fz1m(:,0)*srces(:,l0,1,plus))/vp

      srcm = t1 + t3 + t4

      IF (lscreen) THEN
         PRINT *,'  {F+,sigma+}     = ', srcp/vp
         PRINT *,' -{F-,sigma-}+... = ', srcm
         PRINT *,' -{F-,s-}    = ', t1
         PRINT *,'  {F+,s+}    = ', t2,' {F0+,s+} = ',   t4
         PRINT *,'  {s+,C-1s+} = ', t3
      END IF

      IF (ipmb .eq. 2) THEN
         L11p = 0;         L33p = 0;         L31p = 0;        L13p = 0
         rs11p = 0;        rs33p = 0;        rs13p = 0;       rs31p = 0
         rsd1p = 0;        rsd3p = 0;        crs1p = 0;       crs3p = 0
         stress1 = 0;      bpflux1 = 0
      ELSE
c   See expressions following Eq.(24) in W. I. van Rij, et. al.
         L11p = SUM(fz1p(:,l0)*srces(:,l0,pgrad,plus) +
     1              fz1p(:,l2)*srces(:,l2,pgrad,plus))/vp                !!{f1+, sigma1+}
         rs11p = ABS(g11p/L11p)
         L11p = L11p - g11p                                              !!g11p = {f1+,(sigma1+ -Wf1+)} => 0
         L33p = SUM(fz3p(:,l1)*srces(:,l1,epar,plus) +
     1              fz3p(:,l2)*srces(:,l2,epar,plus))/vp                 !!{f3+, sigma3+}
         rs33p = ABS(g33p/L33p)
         L33p = L33p - g33p
         L13p = SUM(fz1p(:,l1)*srces(:,l1,epar,plus) +
     1              fz1p(:,l2)*srces(:,l2,epar,plus))/vp                 !!{f1+, sigma3+}
         rs13p = ABS(g13p/L13p)
         L13p  = L13p - g13p
         L31p = SUM(fz3p(:,l0)*srces(:,l0,pgrad,plus) +
     1              fz3p(:,l2)*srces(:,l2,pgrad,plus))/vp                !!{f3+, sigma1+}
         rs31p = ABS(g31p/L31p)
         L31p  = L31p - g31p
         stress1 = SUM(fz1p(:,l2)*bstrs)                                 !!l=2, sin component of distribution
         bpflux1 = bpfac*stress1
      ENDIF

      IF (ipmb .eq. 1) THEN
         L11m = 0;     L33m = 0;      L13m = 0;     L31m = 0
         rs11m = 0;    rs33m = 0;     rs13m = 0;    rs31m = 0
         rsd1m  = 0;   rsd3m  = 0;    crs1m  = 0;   crs3m  = 0
      ELSE
         L11m = SUM(fz1m(:,0)*srces(:,l0,pgrad,plus))/vp                 !!{f01+, sigma1+}
         DO mn = l1, l3
            L11m = L11m - SUM(fz1m(:,mn)*srces(:,mn,pgrad,minus))/vp     !!-{f1-, sigma1-}
         END DO
         rs11m = ABS(g11m/L11m)
         L11m  = L11m + g11m - s1cs1                                     !!-{sigma1+, C-1(sigma1+)}

         L33m  =-SUM(fz3m(:,l1)*srces(:,l1,epar,minus))/vp               !!-{f1-, sigma1-}
         rs33m = ABS(g33m/L33m)
         L33m  = L33m + g33m + g33s                                      !!+{Fs, sigma_s}

         L13m  = SUM(fz3m(:,0)*srces(:,l0,pgrad,plus))/vp
         DO mn = l1, l3
            L13m = L13m - SUM(fz3m(:,mn)*srces(:,mn,pgrad,minus))/vp
         END DO
         rs13m = ABS(g13m/L13m)
         L13m  = L13m + g13m

         L31m  =-SUM(fz1m(:,l1)*srces(:,l1,epar,minus))/vp
         rs31m = ABS(g31m/L31m)
         L31m  = L31m + g31m
      ENDIF

      IF (ipmb .ne. 0) THEN
         g31max = 0
         g31min = 0
      ELSE                                                               !!Eq. (26)
         del = .5_dp*SQRT(ABS((L11m - L11p)*(L33m - L33p)))
         av = .25_dp*(L13p + L31p + L13m + L31m)
         g31max = av + del
         g31min = av - del
      ENDIF

c  scale factors relating DIJ in Eq. (36) of W. I. van Rij, et. al.
c  with output from mono-energetic code (these results for LIJ). Note the factor
c  of 0.5 arises due to the conversion from v (velocity) to K (kinetic energy)

      scal11 = 0.5_dp*(b00*(vthermi/wcyclo))**2 * vthermi
      scal33 = 0.5_dp*vthermi
      scal13 = 0.5_dp*(b00*(vthermi/wcyclo)) * vthermi
      
c  output arrays if DKES_rad_dex set.  DKES_rad_dex is set to 0
c  in the main DKES routine, only the STELLOPT routine changes this
c  and allocates these arrays
      IF (ALLOCATED(DKES_L11p) .and. (DKES_rad_dex > 0)) THEN
!         PRINT *,'got here0'
         DKES_L11p(DKES_rad_dex) = L11p
         DKES_L33p(DKES_rad_dex) = L33p
         DKES_L31p(DKES_rad_dex) = L31p
         DKES_L11m(DKES_rad_dex) = L11m
         DKES_L33m(DKES_rad_dex) = L33m
         DKES_L31m(DKES_rad_dex) = L31m
         DKES_scal11(DKES_rad_dex) = scal11
         DKES_scal33(DKES_rad_dex) = scal33
         DKES_scal31(DKES_rad_dex) = scal13
!         PRINT *,'got here1'
      END IF

c  output results summary

      WRITE (ioout, 10) '+', L11p, L33p, L13p, L31p, g31min,
     1                 '-', L11m, L33m, L13m, L31m, g31max,
     2                 scal11, scal33, scal13,
     2                 stress1, bpflux1, g33s,
     3                 'F(1,-)', rsd1m, crs1m, rs11m, rs13m,
     4                 'F(1,+)', rsd1p, crs1p, rs11p, rs13p,
     5                 'F(3,-)', rsd3m, crs3m, rs33m, rs31m,
     6                 'F(3,+)', rsd3p, crs3p, rs33p, rs31p
      CALL FLUSH(ioout)

 10   FORMAT(/1x,'NEOCLASSICAL TRANSPORT MATRIX ELEMENTS: ',
     1     /1x,'DIJ(Eq.36,K=1) = LIJ * [(MKS FACT) * (Ti SCALE)]',
     1      2x,' NOTE: L33(paper) = <B^2>*[L33(Spitzer) - L33]',
     1   //,6x,'PARITY',10x,'L11',12x,'L33',12x,'L13',12x,'L31',7x,
     2      'L13(min/max)',/,1x,131('-'),/,2(9x,a,3x,5(3x,1pe12.4),/),
     3      4x,'MKS FACT.',3(3x,1pe12.4),/,
     4      4x,'Ti SCALE',8x,'TI**1.5',8x,'TI**0.5',8x,'TI**1.0',
     3      //,4x,'PARALLEL',7x,'BAN-PLAT',10x,'L33'/,
     4      6x,'STRESS',9x,'FLUX',9x,'(SPITZER)',/,
     5      1pe12.4,4x,1pe12.4,3x,1pe12.4///,1x,
     5      'EQUATION RESIDUALS (L = KINETIC EQUATION OPERATOR)'//,
     6      3x,' F(I,+/-)',7x,'KIN. EQ.',4x,'PART. CONSERV.',3x,
     7      'RES(FI,LI)',5x,'RES(FJ,LI)',/,
     8      18x,'{L[FI]**2}',20x,'{FI,L[FI]}',5x,'{FJ,L[FI]}',/,
     9      1x,131('-'),/,4(5x,a,2x,4(3x,1pe12.4)/)/)

!
!     WRITE SUMMARY OPT_FILE FOR USE BY OPTIMIZER
!
      WRITE(ioout_opt,'(3(2x,e24.13))') L11p, L33p, L31p
      WRITE(ioout_opt,'(3(2x,e24.13))') L11m, L33m, L31m
      WRITE(ioout_opt,'(3(2x,e24.13))') scal11, scal33, scal13
      CALL FLUSH(ioout_opt)

      !CLOSE(unit=ioout_opt)

c      Generate output for DKES parameter study file            !record file addition
c
c      Find maximim abs(residual):
c
       rsds(1)=rsd1m;rsds(2)=crs1m;rsds(3)=rs11m;rsds(4)=rs13m    !record file addition
       rsds(5)=rsd1p;rsds(6)=crs1p;rsds(7)=rs11p;rsds(8)=rs13p    !record file addition
       rsds(9)=rsd3m;rsds(10)=crs3m;rsds(11)=rs33m;rsds(12)=rs31m !record file addition
       rsds(13)=rsd3p;rsds(14)=crs3p;rsds(15)=rs33p;rsds(16)=rs31p!record file addition
       DO ir=1,16                                                 !record file addition
        rsds(ir) = abs(rsds(ir))                                  !record file addition
       END DO                                                     !record file addition
       rsds_max = maxval(rsds)                                    !record file addition

       WRITE(itab_out,99) cmul1,tb,efield1,tb,weov,tb,wtov,tb,L11m,!record file addition
     >  tb,L11p,tb,L31m,tb,L31p,tb,L33m,tb,L33p,tb,scal11,tb,      !record file addition
     >  scal13,tb,scal33,tb,rsds_max,tb,chip,tb,psip,              !record file addition
     >  tb,btheta,tb,bzeta,tb,vp                                   !record file addition
   99 FORMAT(18(e12.5,a1),e12.5)                                   !record file addition
      !CLOSE(unit=itab_out)                                         !record file addition
      END SUBROUTINE dkes_printout
