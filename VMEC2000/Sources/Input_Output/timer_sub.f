      MODULE timer_sub
      USE stel_kinds, ONLY: dp
      USE parallel_include_module
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, PARAMETER :: tsum = 0, tvac = 1, tread = 2, twout= 3, 
     &                      teqf = 4, tfun = 5, trecon= 6, tfft = 7,
     &                      tffi = 8, tfor = 9, tbcov =10, tres = 11,
     &                      tprec2d = 12, tvac_2d = 13, tfun_2d = 14,
     &                      tfact_2d=15, tio = 16

      INTEGER, PARAMETER :: tsurf=1, tscal=2, tbext=3, tsolver=4,
     &                      tallg=5, tfouri=6, tgreenf=7, tfourp=8,
     &                      tallr=9, tanal=10, tasum=11, tasum2=12,
     &                      tallgv=13, tanar=14
                            
      REAL(dp) :: treadon, treadoff, tfunon, tfunoff,
     &            treconon, treconoff, tffton, tfftoff,
     &            tbcovon, tbcovoff, tvacon, tvacoff,
     &            tforon, tforoff, treson, tresoff,
     &            tprec2don, tprec2doff, twouton, twoutoff,
     &            teqfon, teqfoff, timeon, timeoff,
     &            timer_tsum, timer_tfun, timer_io
      REAL(dp), DIMENSION(0:15) :: timer=0
      REAL(dp), DIMENSION(15) :: timer_vac=0

      CONTAINS
      
      SUBROUTINE write_times (nthreed, lscreen, lfreeb, lrecon, lprec2d)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: nthreed
      LOGICAL, INTENT(in) :: lscreen, lfreeb, lrecon, lprec2d
      INTEGER             :: i, nform
      CHARACTER(LEN=*), DIMENSION(0:16), PARAMETER :: form =
     &    (/ 'TOTAL COMPUTATIONAL TIME (SEC) ',
     &       '   FREE BOUNDARY (VACUUM)      ',
     &       '   READ IN DATA                ',
     &       '   WRITE OUT DATA TO WOUT      ',
     &       '   EQFORCE                     ',
     &       'TIME IN FUNCT3D                ',
     &       '   PROFILE RECONSTRUCTION      ',
     &       '   FOURIER TRANSFORM           ',
     &       '   INVERSE FOURIER TRANSFORM   ',
     &       '   FORCES AND SYMMETRIZE       ',
     &       '   BCOVAR FIELDS               ',
     &       '   RESIDUE                     ',
     &       'TIME IN PRECON2D SETUP         ',
     &       '   VACUUM ONLY                 ',
     &       '   FUNCT3D                     ',
     &       '   FORWARD SOLVE (FACTOR BLKS) ',
     &       'TIME TO INPUT/OUTPUT           '
     &     /)

      timer_tsum = timer(tsum) + timer(twout) + timer(teqf) 
      timer_tfun = timer(tfun)
      timer_io   = timer(tread) + timer(twout)

      DO i = 1,2
         IF (i .eq. 1) nform = 6
         IF (i .eq. 2) nform = nthreed
         IF (.not.lscreen .and. i.eq.1) CYCLE
         WRITE (nform, 20)
     &         form(tsum) ,timer_tsum,   form(tio), timer_io,
     &         form(tread),timer(tread), form(twout),timer(twout),
     &         form(tfun) , timer(tfun), form(tbcov) ,timer(tbcov),
     &         form(tfft) , timer(tfft), form(tffi)  ,timer(tffi),
     &         form(tfor) , timer(tfor), form(tres), timer(tres),
     &         form(teqf) ,timer(teqf)
         IF (lrecon) WRITE (nform, 20) form(trecon),timer(trecon)
         IF (lfreeb) THEN
            WRITE (nform, 20) form(tvac) ,timer(tvac)
            WRITE (nform, 24) timer_vac(tsurf), timer_vac(tbext), 
     &      timer_vac(tscal), timer_vac(tanal), timer_vac(tasum),
     &      timer_vac(tasum2), timer_vac(tanar), timer_vac(tgreenf),
     &      timer_vac(tfourp), timer_vac(tallr), timer_vac(tallg),
     &      timer_vac(tfouri), timer_vac(tallgv),timer_vac(tsolver)
         END IF
         IF (lprec2d) THEN
	     WRITE (nform, 20) form(tprec2d),  timer(tprec2d),
     &                    form(tfun_2d),  timer(tfun_2d),
     &                    form(tfact_2d), timer(tprec2d)-timer(tfun_2d)
            IF (lfreeb) WRITE (nform, 20) form(tvac_2d), timer(tvac_2d)
         END IF
      END DO

   20 FORMAT(a35,f12.2)
   24 FORMAT(  10x, 'VACUUM SURFACE    ',7x,f12.2,
     &       /,10x, 'VACUUM BEXTERN    ',7x,f12.2,
     &       /,10x, 'VACUUM SCALPOT    ',7x,f12.2,
     &       /,10x, '   ANALYT         ',7x,f12.2,
     &       /,10x, '      ASUM        ',7x,f12.2,
     &       /,10x, '      ASUM2       ',7x,f12.2,
     &       /,10x, '      ALLREDUCE   ',7x,f12.2,
     &       /,10x, '   GREENF         ',7x,f12.2,
     &       /,10x, '   FOURP          ',7x,f12.2,
     &       /,10x, '   ALLREDUCE      ',7x,f12.2,
     &       /,10x, '   ALLGATHER      ',7x,f12.2,
     &       /,10x, '   FOURI          ',7x,f12.2,
     &       /,10x, 'VACUUM ALLGATHER  ',7x,f12.2,
     &       /,10x, 'VACUUM SOLVER     ',7x,f12.2)

      END SUBROUTINE write_times
    
      END MODULE timer_sub
