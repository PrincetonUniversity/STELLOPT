      SUBROUTINE getbpsq(psixz,nxd,nzd,xgrid,dx,dz,nx,nz,bpsq)
      USE precision
      IMPLICIT NONE
      INTEGER :: nxd,nzd,nx,nz
      REAL(rprec) :: psixz(nxd,nzd)
      REAL(rprec) :: bpsq(nxd,nzd)
      REAL(rprec) :: xgrid(nxd)
      REAL(rprec) :: dx,dz
      REAL(rprec) :: dpsixsq
      INTEGER :: i,j
      DO i=2,nx-1
         DO j=2,nz-1
            bpsq(i,j)=(((psixz(i+1,j)-psixz(i-1,j))/(2.*dx))**2+
     $          ((psixz(i,j+1)-psixz(i,j-1))/(2.*dz))**2)/
     $           xgrid(i)**2
         END DO
      END DO
      i=1
      DO j=2,nz-1
         bpsq(i,j)=(((-3.*psixz(i,j)+4.*psixz(i+1,j)
     $        -psixz(i+2,j))/(2.*dx))**2+
     $        ((psixz(i,j+1)-psixz(i,j-1))/(2.*dz))**2)/
     $        xgrid(i)**2
      END DO
      i=nx
      DO j=2,nz-1
         bpsq(i,j)=(((3.*psixz(i,j)-4.*psixz(i-1,j)
     $        +psixz(i-2,j))/(2.*dx))**2+
     $        ((psixz(i,j+1)-psixz(i,j-1))/(2.*dz))**2)/
     $        xgrid(i)**2
      END DO
      j=1
      DO i=1,nx
         IF(i.eq.1) THEN
            dpsixsq=((-3.*psixz(i,j)+4.*psixz(i+1,j)
     $        -psixz(i+2,j))/(2.*dx))**2
         ELSE IF(i.eq.nx) THEN
            dpsixsq=((3.*psixz(i,j)-4.*psixz(i-1,j)
     $        +psixz(i-2,j))/(2.*dx))**2
         ELSE
            dpsixsq=((psixz(i+1,j)-psixz(i-1,j))/(2.*dx))**2
         ENDIF
         bpsq(i,j)=(dpsixsq+
     $        ((-3.*psixz(i,j)+4.*psixz(i,j+1)
     $        -psixz(i,j+2))/(2.*dz))**2)/
     $        xgrid(i)**2
      END DO
      j=nz
      DO i=1,nx
         IF(i.eq.1) THEN
            dpsixsq=((-3.*psixz(i,j)+4.*psixz(i+1,j)
     $        -psixz(i+2,j))/(2.*dx))**2
         ELSE IF(i.eq.nx) THEN
            dpsixsq=((3.*psixz(i,j)-4.*psixz(i-1,j)
     $        +psixz(i-2,j))/(2.*dx))**2
         ELSE
            dpsixsq=((psixz(i+1,j)-psixz(i-1,j))/(2.*dx))**2
         ENDIF
         bpsq(i,j)=(dpsixsq+
     $        ((3.*psixz(i,j)-4.*psixz(i,j-1)
     $        +psixz(i,j-2))/(2.*dz))**2)/
     $        xgrid(i)**2
      END DO
      RETURN
      END
