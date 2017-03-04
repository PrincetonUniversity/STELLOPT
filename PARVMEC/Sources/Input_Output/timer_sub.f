      MODULE timer_sub
      USE stel_kinds, ONLY: rprec
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
     3                      tprec2d = 12

!SPH CHANGE
      INTEGER, PARAMETER :: tsurf=1, tscal=2, tbext=3, tsolver=4,
     1                      tallg=5, tfouri=6, tgreenf=7, tfourp=8,
     2                      tallr=9, tanal=10, tasum=11, tasum2=12, 
     3                      tallgv=13, tanar=14, tbcast_vac=15
                            
      REAL(rprec) :: tvacon, tvacoff, tfunon, tfunoff,
     1               twouton, twoutoff, teqfon, teqfoff,
     2               treadon, treadoff, timeon, timeoff,
     3               treconon, treconoff, tffton, tfftoff,
     4               tbcovon, tbcovoff, tforon, tforoff,
     5               treson, tresoff, tprec2don, tprec2doff,
     6               timer_tsum, timer_tfun
      REAL(rprec), DIMENSION(0:12) :: timer
!SPH CHANGE
      REAL(rprec), DIMENSION(15) :: timer_vac=0

      CONTAINS
      
      SUBROUTINE write_times (nthreed, lscreen, lfreeb, lrecon, lprec2d)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: nthreed
      LOGICAL, INTENT(in) :: lscreen, lfreeb, lrecon, lprec2d
      INTEGER             :: i, nform
      CHARACTER(LEN=*), DIMENSION(0:12), PARAMETER :: form =
     1    (/ 'TOTAL COMPUTATIONAL TIME       ',
     2       'TIME IN VACUUM LOOP            ',
     3       'TIME TO READ IN DATA           ',
     4       'TIME TO WRITE DATA TO WOUT     ',
     5       'TIME IN EQFORCE                ',
     6       'TIME (REMAINDER) IN FUNCT3D    ',
     7       'TIME IN PROFILE RECONSTRUCTION ',
     8       'TIME IN FOURIER TRANSFORM      ',
     9       'TIME IN INVERSE FOURIER XFORM  ',
     A       'TIME IN FORCES + SYMFORCES     ',
     B       'TIME IN BCOVAR                 ',
     C       'TIME IN RESIDUE                ',
     D       'TIME IN PRECON2D SETUP         '
     Z     /)

      timer_tsum = timer(tsum) + timer(twout) + timer(teqf) 

      timer_tfun = timer(tfun) - timer(tvac) - timer(trecon)
     1           - timer(tfft) - timer(tffi) - timer(tfor) 
     2           - timer(tbcov) - timer(tres)

      DO i = 1,2
         IF (i .eq. 1) nform = 6
         IF (i .eq. 2) nform = nthreed
         IF (.not.lscreen .and. i.eq.1) CYCLE
         IF (lfreeb .and. lrecon) THEN
            WRITE (nform, 20)
     1            form(tsum) ,timer_tsum  , form(tread) ,timer(tread), 
     2            form(twout),timer(twout), form(teqf)  ,timer(teqf),
     3            form(tvac) ,timer(tvac) , form(trecon),timer(trecon),
     4            form(tfft), timer(tfft),  form(tffi)  ,timer(tffi),
     5            form(tfor), timer(tfor),  form(tbcov) ,timer(tbcov),
     6            form(tres), timer(tres)
         ELSE IF (lfreeb) THEN
            WRITE (nform, 20)
     1            form(tsum),timer_tsum,form(tread),timer(tread),
     2            form(twout),timer(twout),form(teqf),timer(teqf), 
     3            form(tvac),timer(tvac),
     4            form(tfft), timer(tfft), form(tffi), timer(tffi),
     5            form(tfor), timer(tfor), form(tbcov),timer(tbcov),
     6            form(tres), timer(tres)
         ELSE IF (lrecon) THEN
            WRITE (nform, 20) 
     1            form(tsum), timer_tsum, form(tread),timer(tread),
     2            form(twout),timer(twout),form(teqf), timer(teqf),
     3            form(trecon),timer(trecon),
     4            form(tfft), timer(tfft), form(tffi), timer(tffi),
     5            form(tfor), timer(tfor), form(tbcov),timer(tbcov),
     6            form(tres), timer(tres)
         ELSE
            WRITE (nform, 20) 
     1            form(tsum),timer_tsum,form(tread),timer(tread),
     2            form(twout),timer(twout),form(teqf),timer(teqf), 
     3            form(tfft), timer(tfft), form(tffi), timer(tffi),
     4            form(tfor), timer(tfor), form(tbcov),timer(tbcov),
     5            form(tres), timer(tres) 
         END IF
         IF (lprec2d)  WRITE (nform, 20) form(tprec2d), timer(tprec2d)
         WRITE (nform, 20) form(tfun), timer_tfun
!SPH CHANGE
         IF (lfreeb) THEN
            WRITE (nform, 24)
     1      timer_vac(tbcast_vac), timer_vac(tsurf), timer_vac(tbext), 
     2      timer_vac(tscal), timer_vac(tanal), timer_vac(tasum), 
     3      timer_vac(tasum2), timer_vac(tanar), timer_vac(tgreenf), 
     4      timer_vac(tfourp), timer_vac(tallr), timer_vac(tallg),
     5      timer_vac(tfouri), timer_vac(tallgv),timer_vac(tsolver)
         ENDIF
      END DO

   20 FORMAT(a35,f12.2,' SECONDS')
   24 FORMAT(/,4x, 'FUNCT3D BCAST     ',f12.2,
     1       /,4x, 'VACUUM SURFACE    ',f12.2,
     1       /,4x, 'VACUUM BEXTERN    ',f12.2,
     2       /,4x, 'VACUUM SCALPOT    ',f12.2,
     3       /,4x, '   ANALYT         ',f12.2,
     3       /,4x, '      ASUM        ',f12.2,
     3       /,4x, '      ASUM2       ',f12.2,
     3       /,4x, '      ALLREDUCE   ',f12.2,
     3       /,4x, '   GREENF         ',f12.2,
     3       /,4x, '   FOURP          ',f12.2,
     3       /,4x, '   ALLREDUCE      ',f12.2,
     3       /,4x, '   ALLGATHER      ',f12.2,
     3       /,4x, '   FOURI          ',f12.2,
     3       /,4x, 'VACUUM ALLGATHER  ',f12.2,
     4       /,4x, 'VACUUM SOLVER     ',f12.2)

      END SUBROUTINE write_times
    
      END MODULE timer_sub
