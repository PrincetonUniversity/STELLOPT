      SUBROUTINE write_gade_nml(iunit)
      USE stel_kinds
      USE gade_mod
      USE write_array_generic
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: iunit
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i
C-----------------------------------------------
      WRITE (iunit,'(a)')'&GA_DE'
      WRITE (iunit,200) 'NPOPSIZ = ', npopsiz
      WRITE (iunit,200) 'NGEN    = ', ngen
      WRITE (iunit,200) 'IDUM    = ', idum
      WRITE (iunit,200) 'IBOUND  = ', ibound
      WRITE (iunit,200) 'NOWRITE = ', nowrite
      WRITE (iunit,200) 'MICROGA = ', microga
      WRITE (iunit,200) 'ISKIP   = ', iskip
      WRITE (iunit,200) 'IEND    = ', iend
      WRITE (iunit,200) 'NCHILD  = ', nchild
      WRITE (iunit,200) 'ITOURNY = ', itourny
      WRITE (iunit,200) 'IELITE  = ', ielite
      WRITE (iunit,200) 'ICREEP  = ', icreep
      WRITE (iunit,200) 'IUNIFRM = ', iunifrm
      WRITE (iunit,200) 'INICHE  = ', iniche
      WRITE (iunit,200) 'STRATEGY= ', strategy
      WRITE (iunit,200) 'CR_STRATEGY = ', cr_strategy
      WRITE (iunit,200) 'OUT_ITER= ', out_iter
      WRITE (iunit,200) 'UNIQUE_IND = ', unique_ind
      WRITE (iunit,210) 'PCROSS  = ', pcross
      WRITE (iunit,210) 'F_CROSS = ', f_cross
      WRITE (iunit,210) 'PMUTATE = ', pmutate
      WRITE (iunit,210) 'PCREEP  = ', pcreep
      WRITE (iunit,101) 'SAVE_SPACE = ', save_space
      WRITE (iunit,100) 'NICHFLG = '
      WRITE (iunit,110) (nichflg(i), i=1,nparmax)
      WRITE (iunit,100) 'NPOSIBL = '
      WRITE (iunit,110) (nposibl(i), i=1,nparmax)
      CALL write_array (iunit, 'PARMIN', parmin, nparmax)
      CALL write_array (iunit, 'PARMAX', parmax, nparmax)
      WRITE (iunit,'(a)') '/'

 100  FORMAT (2x,a)
 101  FORMAT (2x,a,l2)
 110  FORMAT (2x, 16i5)
 200  FORMAT (2x,a,i5)
 210  FORMAT (2x,a,1pe21.14)

      END SUBROUTINE write_gade_nml
