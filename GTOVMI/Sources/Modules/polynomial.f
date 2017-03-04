      MODULE polynomial
      USE precision
      CONTAINS
        FUNCTION polyval(a,n,x,m)
         IMPLICIT NONE
         INTEGER :: n, m, nt, i, j
         REAL(rprec) a(0:n-1),x(*)
         REAL(rprec), POINTER :: polyval(:)
         nt=n-1
         ALLOCATE(polyval(m))
         DO i=1,m
           polyval(i)=a(0)
           DO j=1,nt
              IF(x(i).ne.0.)polyval(i)=polyval(i)+a(j)*x(i)**j
           ENDDO
         ENDDO
         RETURN
        END FUNCTION polyval
        FUNCTION polyder(a,n,x,m)
         IMPLICIT NONE
         INTEGER :: n, m, nt, i, j
         REAL(rprec) a(0:n-1),x(*)
         REAL(rprec), POINTER :: polyder(:)
         nt=n-1
         ALLOCATE(polyder(m))
         DO i=1,m
           polyder(i)=a(1)
           DO j=2,nt
              IF(x(i).ne.0.)polyder(i)=polyder(i)+a(j)*j*x(i)**(j-1)
           ENDDO
         ENDDO
         RETURN
        END FUNCTION polyder
        FUNCTION polyint(a,n,x,m)
         IMPLICIT NONE
         INTEGER :: n, m, nt, i, j
         REAL(rprec) a(0:n-1),x(*)
         REAL(rprec), POINTER :: polyint(:)
         nt=n-1
         ALLOCATE(polyint(m))
         polyint=0.
         DO i=1,m
           polyint(i)=a(0)*x(i)
           DO j=1,nt
              IF(x(i).ne.0.)polyint(i)=polyint(i)+a(j)/(j+1)*x(i)**(j+1)
           ENDDO
         ENDDO
         RETURN
        END FUNCTION polyint

      END MODULE polynomial
