      SUBROUTINE printout(i0, delt0, w0, lscreen)
      USE vmec_main
      USE realspace
      USE xstuff
#ifdef _HBANGLE
      USE angle_constraints, ONLY: getrz
#endif
      USE parallel_include_module
      USE parallel_vmec_module, ONLY: CopyLastNtype
      USE vmec_params, ONLY: ntmax
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: i0
      REAL(dp) :: delt0, w0
      LOGICAL :: lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      CHARACTER(LEN=*), PARAMETER :: 
#ifdef _HBANGLE
     1   iter_line = "  ITER    FSQRHO   FSQ(m=0)    FSQL   ",
     1   iter_lines= "  ITER  FSQRHO   FSQ(m=0)    FSQL   ",
     2   fsq_line  = "  fsqrho   fsq(m=0)    fsql      DELT    ",
     2   fsq_lines = "  fsqrho   fsq(m=0)    fsql      DELT    ",
     4   raxis_line = " RAX(v=0) ",
#else
     1   iter_line = "  ITER    FSQR      FSQZ      FSQL   ",
     2   fsq_line  = "   fsqr      fsqz      fsql      DELT    ",
     3   iter_lines = iter_line, fsq_lines = fsq_line,
     4   raxis_line = "RAX(v=0) ",
#endif
     3   delt_line = "    DELT  ",        !J.Geiger 20101118
     5   zaxis_line = "  ZAX(v=0)      "
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(dp) :: betav, w, avm, den, tbroadon, tbroadoff
      REAL(dp), ALLOCATABLE :: bcastbuf(:)
      CHARACTER(len=LEN(iter_line) + LEN(fsq_line) +
     1          LEN(raxis_line) + LEN(zaxis_line)) :: 
     2          print_line
      INTEGER :: i, j, k, l, lk
C-----------------------------------------------
      IF(grank .GE. nranks) RETURN

#ifdef _ANIMEC
      betav = (2*wper + wpar)/(3*wb)
#else
      betav = wp/wb
#endif
      w = w0*twopi*twopi
      den = zero
      specw(1) = one

      IF(PARVMEC) THEN 
        CALL CopyLastNtype(pxstore, pgc)
      ELSE
        gc = xstore
      END IF
#ifdef _HBANGLE
      CALL getrz(gc)
#endif 

      IF (PARVMEC) THEN 
        CALL spectrum_par (pgc(:irzloff), pgc(1+irzloff:2*irzloff))
        CALL Gather1XArray(vp)
        CALL Gather1XArray(specw)
      ELSE
        CALL spectrum (gc(:irzloff), gc(1+irzloff:2*irzloff))
      END IF

      den = SUM(vp(2:ns))
      avm = DOT_PRODUCT(vp(2:ns),specw(2:ns)+specw(1:ns-1))
      avm = 0.5_dp*avm/den
      IF (ivac .GE. 1) THEN
        IF (PARVMEC) THEN
!SPH CHANGE (MOVE OUT OF FUNCT3D)
#if defined(MPI_OPT)
          ACTIVE1: IF (lactive) THEN
          CALL second0(tbroadon)
          ALLOCATE(bcastbuf(3*nznt+1))
          bcastbuf(1:nznt) = dbsq
          bcastbuf(nznt+1:2*nznt) = bsqsav(:,3)
          bcastbuf(2*nznt+1:3*nznt) = rbsq    !NEED THIS WHEN INTERPOLATING MESHES
          bcastbuf(3*nznt+1) = fedge

          CALL MPI_Bcast(bcastbuf,SIZE(bcastbuf),MPI_REAL8,
     1                   nranks-1,NS_COMM,MPI_ERR)

          dbsq = bcastbuf(1:nznt)
          bsqsav(:,3) = bcastbuf(nznt+1:2*nznt)
          rbsq = bcastbuf(2*nznt+1:3*nznt)
          fedge = bcastbuf(3*nznt+1)
          DEALLOCATE(bcastbuf)

          CALL second0(tbroadoff)
          broadcast_time = broadcast_time + (tbroadoff - tbroadon)
          den = SUM(bsqsav(:nznt,3)*pwint(:,2))
          IF (den .NE. zero) delbsq =
     1     SUM(dbsq(:nznt)*pwint(:,2))/den
          END IF ACTIVE1
#endif
        ELSE
          delbsq =
     1     SUM(dbsq(:nznt)*wint(2:nrzt:ns))/
     2     SUM(bsqsav(:nznt,3)*wint(2:nrzt:ns))
        END IF
      END IF

      IF (i0.EQ.1 .AND. lfreeb) THEN
         print_line = iter_lines // " " // raxis_line 
         IF (lasym) print_line = TRIM(print_line) // " " // zaxis_line
         IF (lscreen.AND.grank.EQ.0) 
     1        PRINT 20, TRIM(print_line)//delt_line  !J Geiger 20101118
         print_line = iter_line // fsq_line // raxis_line
         IF (lasym) print_line = TRIM(print_line) // " " // zaxis_line
         IF (imatch_phiedge .eq. 1) THEN
           IF(grank.EQ.0) WRITE (nthreed, 15) TRIM(print_line)
         ELSE
           IF(grank.EQ.0) WRITE (nthreed, 16) TRIM(print_line)
         ENDIF
      ELSE IF (i0.eq.1 .and. .not.lfreeb) THEN
         print_line = raxis_line 
         IF (lasym) print_line = raxis_line // zaxis_line
         IF (lscreen.AND.grank.EQ.0) 
     1      PRINT 30, iter_lines, TRIM(print_line)//delt_line !J Geiger 2010118
         print_line = iter_line // fsq_line // raxis_line // "     "
         IF (lasym) print_line = iter_line // fsq_line // raxis_line
     1                        // zaxis_line
         IF (grank .EQ. 0) WRITE (nthreed, 25) TRIM(print_line)
      ENDIF
   15 FORMAT(/,a,6x,'WMHD      BETA      <M>   DEL-BSQ   FEDGE',/)
   16 FORMAT(/,a,6x,'WMHD      BETA     PHIEDGE  DEL-BSQ    FEDGE',/)
   20 FORMAT(/,a,6x,'WMHD      DEL-BSQ',/)
   25 FORMAT(/,a,6x,'WMHD      BETA      <M>        ',/)
   30 FORMAT(/,a,1x,a,5x,'WMHD',/)

      IF (.not. lasym) THEN
         IF (.not.lfreeb) THEN
            IF (lscreen.AND.grank.EQ.0) 
     1        PRINT 45, i0, fsqr, fsqz, fsql, r00, delt0, w !J Geiger 20101118
            IF(grank.EQ.0) WRITE (nthreed, 40) i0, fsqr, fsqz, fsql,
     1      fsqr1, fsqz1, fsql1, delt0, r00, w, betav, avm
            RETURN
         ENDIF
         IF (lscreen.AND.grank.EQ.0) 
     1        PRINT 50, i0, fsqr, fsqz, fsql, r00, delt0, w,
     2                          delbsq !J Geiger 20101118
         IF (imatch_phiedge .eq. 1) THEN
            IF(grank.EQ.0) WRITE (nthreed, 40) i0, fsqr, fsqz, fsql,
     1      fsqr1, fsqz1,  fsql1, delt0, r00, w, betav, avm, delbsq,
     2      fedge
         ELSE
            IF(grank.EQ.0) WRITE (nthreed, 42) i0, fsqr, fsqz, fsql,
     1      fsqr1, fsqz1, fsql1, delt0, r00, w, betav,
     2      ABS(phiedge), delbsq, fedge
         ENDIF
      
      ELSE
         IF (.not.lfreeb) THEN
            IF (lscreen.AND.grank.EQ.0) 
     1        PRINT 65, i0, fsqr, fsqz, fsql, r00, z00, !J Geiger 20101118
     2                             delt0, w !J Geiger 20101118
            IF(grank.EQ.0) WRITE (nthreed, 60) i0, fsqr, fsqz, fsql,
     1       fsqr1, fsqz1, fsql1, delt0, r00, z00, w, betav, avm
         RETURN
         ENDIF
         IF (lscreen.AND.grank.EQ.0) 
     1        PRINT 70, i0, fsqr, fsqz, fsql, r00, z00,
     2                          delt0, w, delbsq !J Geiger 20101118
         IF (imatch_phiedge .eq. 1) THEN
            IF(grank.EQ.0) WRITE (nthreed, 60) i0, fsqr, fsqz, fsql,
     1       fsqr1, fsqz1, fsql1, delt0, r00, z00, w, betav, avm,
     2        delbsq, fedge
         ELSE
            IF(grank.EQ.0) WRITE (nthreed, 60) i0, fsqr, fsqz, fsql,
     1      fsqr1, fsqz1, fsql1, delt0, r00, z00, w, betav, 
     2      ABS(phiedge), delbsq, fedge
         ENDIF
      END IF

   40 FORMAT(i6,1x,1p,7e10.2,e11.3,e12.4,e11.3,0p,f7.3,1p,2e9.2)
   42 FORMAT(i5,1p,7e10.2,e11.3,e12.4,2e11.3,0p,f7.3,1p,e9.2)
   45 FORMAT(i5,1p,3e10.2,e11.3,e10.2,e12.4)
   50 FORMAT(i5,1p,3e10.2,e11.3,e10.2,e12.4,e11.3)
   60 FORMAT(i6,1x,1p,7e10.2,2e11.3,e12.4,e11.3,0p,f7.3,1p,2e9.2)
   65 FORMAT(i5,1p,3e10.2,2e11.3,e10.2,e12.4)
   70 FORMAT(i5,1p,3e10.2,2e11.3,e10.2,e12.4,e11.3)

      END SUBROUTINE printout
