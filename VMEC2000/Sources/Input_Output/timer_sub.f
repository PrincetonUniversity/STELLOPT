      MODULE timer_sub
      USE stel_kinds, ONLY: dp
#if defined(SKS)
      USE parallel_include_module
#endif      
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, PARAMETER :: tsum = 0, tvac = 1, tread = 2, twout= 3, 
     1                      teqf = 4, tfun = 5, trecon= 6, tfft = 7,
     2                      tffi = 8, tfor = 9, tbcov =10, tres = 11,
     3                      tprec2d = 12, tvac_2d = 13, tfun_2d = 14,
     4                      tfact_2d=15, tio = 16

      INTEGER, PARAMETER :: tsurf=1, tscal=2, tbext=3, tsolver=4,
     1                      tallg=5, tfouri=6, tgreenf=7, tfourp=8,
     2                      tallr=9, tanal=10, tasum=11, tasum2=12, 
     3                      tallgv=13, tanar=14
                            
      REAL(dp) :: treadon, treadoff, tfunon, tfunoff,
     1            treconon, treconoff, tffton, tfftoff,
     2            tbcovon, tbcovoff, tvacon, tvacoff,
     3            tforon, tforoff, treson, tresoff, 
     4            tprec2don, tprec2doff, twouton, twoutoff, 
     5            teqfon, teqfoff, timeon, timeoff,
     6            timer_tsum, timer_tfun, timer_io
      REAL(dp), DIMENSION(0:15) :: timer=0
      REAL(dp), DIMENSION(15) :: timer_vac=0

      CONTAINS
      
      SUBROUTINE write_times (nthreed, lscreen, lfreeb, lrecon, lprec2d)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: nthreed
      LOGICAL, INTENT(in) :: lscreen, lfreeb, lrecon, lprec2d
      INTEGER             :: i, nform
      CHARACTER(LEN=*), DIMENSION(0:16), PARAMETER :: form =
     1    (/ 'TOTAL COMPUTATIONAL TIME (SEC) ',
     2       '   FREE BOUNDARY (VACUUM)      ',
     3       '   READ IN DATA                ',
     4       '   WRITE OUT DATA TO WOUT      ',
     5       '   EQFORCE                     ',
     6       'TIME IN FUNCT3D                ',
     7       '   PROFILE RECONSTRUCTION      ',
     8       '   FOURIER TRANSFORM           ',
     9       '   INVERSE FOURIER TRANSFORM   ',
     A       '   FORCES AND SYMMETRIZE       ',
     B       '   BCOVAR FIELDS               ',
     C       '   RESIDUE                     ',
     D       'TIME IN PRECON2D SETUP         ',
     E       '   VACUUM ONLY                 ',
     F       '   FUNCT3D                     ',
     G       '   FORWARD SOLVE (FACTOR BLKS) ',
     H       'TIME TO INPUT/OUTPUT           '
     Z     /)

      timer_tsum = timer(tsum) + timer(twout) + timer(teqf) 
      timer_tfun = timer(tfun)
      timer_io   = timer(tread) + timer(twout)

      DO i = 1,2
         IF (i .eq. 1) nform = 6
         IF (i .eq. 2) nform = nthreed
         IF (.not.lscreen .and. i.eq.1) CYCLE
         WRITE (nform, 20)
     1         form(tsum) ,timer_tsum,   form(tio), timer_io,
     2         form(tread),timer(tread), form(twout),timer(twout), 
     3         form(tfun) , timer(tfun), form(tbcov) ,timer(tbcov),
     4         form(tfft) , timer(tfft), form(tffi)  ,timer(tffi),
     5         form(tfor) , timer(tfor), form(tres), timer(tres),
     6         form(teqf) ,timer(teqf)
         IF (lrecon) WRITE (nform, 20) form(trecon),timer(trecon)
         IF (lfreeb) THEN
            WRITE (nform, 20) form(tvac) ,timer(tvac)
            WRITE (nform, 24) timer_vac(tsurf), timer_vac(tbext), 
     1      timer_vac(tscal), timer_vac(tanal), timer_vac(tasum), 
     2      timer_vac(tasum2), timer_vac(tanar), timer_vac(tgreenf), 
     3      timer_vac(tfourp), timer_vac(tallr), timer_vac(tallg),
     4      timer_vac(tfouri), timer_vac(tallgv),timer_vac(tsolver)
         END IF
         IF (lprec2d) THEN
	    WRITE (nform, 20) form(tprec2d),  timer(tprec2d),
     1                    form(tfun_2d),  timer(tfun_2d),
     2                    form(tfact_2d), timer(tprec2d)-timer(tfun_2d)
            IF (lfreeb) WRITE (nform, 20) form(tvac_2d), timer(tvac_2d)
         END IF
      END DO

   20 FORMAT(a35,f12.2)
   24 FORMAT(  10x, 'VACUUM SURFACE    ',7x,f12.2,
     1       /,10x, 'VACUUM BEXTERN    ',7x,f12.2,
     2       /,10x, 'VACUUM SCALPOT    ',7x,f12.2,
     3       /,10x, '   ANALYT         ',7x,f12.2,
     3       /,10x, '      ASUM        ',7x,f12.2,
     3       /,10x, '      ASUM2       ',7x,f12.2,
     3       /,10x, '      ALLREDUCE   ',7x,f12.2,
     3       /,10x, '   GREENF         ',7x,f12.2,
     3       /,10x, '   FOURP          ',7x,f12.2,
     3       /,10x, '   ALLREDUCE      ',7x,f12.2,
     3       /,10x, '   ALLGATHER      ',7x,f12.2,
     3       /,10x, '   FOURI          ',7x,f12.2,
     3       /,10x, 'VACUUM ALLGATHER  ',7x,f12.2,
     4       /,10x, 'VACUUM SOLVER     ',7x,f12.2)

      END SUBROUTINE write_times
    
      END MODULE timer_sub

