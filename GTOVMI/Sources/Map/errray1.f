      SUBROUTINE errray1(errcode,psi,thet,x,y,err,tol,tola)
      USE precision
      IMPLICIT NONE
      INTEGER :: errcode
      REAL(rprec) :: psi,thet,x,y,err,tol,tola
c.......................................................................
      IF (errcode .eq. 1)then
         WRITE (6,'(a)')'rayatt:  sliding window search out of range'
         WRITE (6,100)'  psi = ',psi,'  thet = ',thet
         WRITE (6,100)'  x   = ',x,  '  y    = ',y
         RETURN
      ENDIF
c
      IF (errcode .eq. 2)then
         WRITE (6,'(a)')'rayatt:  sliding window search fails (10)'
         WRITE (6,100)'  psi = ',psi,'  thet = ',thet
         WRITE (6,100)'  x   = ',x,  '  y    = ',y
         RETURN
      ENDIF
c
      IF (errcode .eq. 3)then
         WRITE (6,'(a)')'rayatt:  Newton-Raphson INTerations fails (20)'
         WRITE (6,100)'  psi = ',psi,'  thet = ',thet
         WRITE (6,100)'  err = ',err,'  tol  = ',tol
         IF (ABS(err) .gt. ABS(tola))then
            WRITE (6,'(a)')'rayatt: err > max.tol'
            RETURN
         ELSE
            WRITE (6,'(a)')'rayatt: tol < err < max. tol (warning)'
            RETURN
         ENDIF
      ENDIF
c
      WRITE(6,*) 'error in errray1--errray999'
      STOP
 100  FORMAT(a,1pe12.4,a,1pe12.4)
      END
