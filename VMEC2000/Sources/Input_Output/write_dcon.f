      SUBROUTINE write_dcon (rzl_array)
      USE vmec_main, fpsi=>bvcof
      USE vmec_params, ONLY: ntmax, rcc, rsc, zsc, zcc, mscale, nscale
      USE realspace
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax),
     1   TARGET, INTENT(in) :: rzl_array
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: istat, m
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: rmncc, rmnsc, zmnsc,    &
     &   zmncc, lmnsc, lmncc
      REAL(rprec) :: t1(0:mpol1)
      CHARACTER*(256) :: dcon_file
C-----------------------------------------------
      t1 = nscale(0)*mscale(0:mpol1)

      ALLOCATE (rmncc(ns,0:0,0:mpol1), zmnsc(ns,0:0,0:mpol1),           &
     &          lmnsc(ns,0:0,0:mpol1), stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error in write_dcon'
      DO m=0,mpol1
         rmncc(:,0,m) = rzl_array(:,0,m,rcc)*t1(m)               !!COS(mu) COS(nv)
         zmnsc(:,0,m) = rzl_array(:,0,m,zsc+ntmax)*t1(m)         !!SIN(mu) COS(nv)
         lmnsc(:,0,m) = rzl_array(:,0,m,zsc+2*ntmax)*t1(m)       !!SIN(mu) COS(nv)
      END DO

      IF (lasym) THEN
      ALLOCATE (rmnsc(ns,0:0,0:mpol1), zmncc(ns,0:0,0:mpol1),           &
     &          lmncc(ns,0:0,0:mpol1), stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error in write_dcon'
      DO m=0,mpol1
         rmnsc(:,0,m) = rzl_array(:,0,m,rsc)*t1(m)               !!SIN(mu) COS(nv)
         zmncc(:,0,m) = rzl_array(:,0,m,zcc+ntmax)*t1(m)         !!COS(mu) COS(nv)
         lmncc(:,0,m) = rzl_array(:,0,m,zcc+2*ntmax)*t1(m)       !!COS(mu) COS(nv)
      END DO
      ENDIF

!     HERE, FULL(i) = (i-1)*hs, i=1,ns (hs=1/(ns-1))
!            HALF(i) = (i-1.5)*hs, i=2,ns

      dcon_file = "dcon_" // TRIM(input_extension) // ".txt"
      OPEN (unit=51,FILE=dcon_file,FORM='FORMATTED',iostat=istat)
      IF (istat .ne. 0) STOP 'Error writing dcon input file'

      IF (mnmax .ne. mpol) STOP 'THIS IS NOT AXISYMMETRIC!'

      WRITE (51, *) ns                             !Number of flux surfaces
      WRITE (51, *) mpol                           !Number of poloidal modes, m=[0:mpol-1]
      WRITE (51, *) lasym                          !Up-down sym:=F
      WRITE (51, *) rmncc(1:ns,0,0:mpol1)          !r = sum [rmnc * cos (mu)], full mesh
      WRITE (51, *) zmnsc(1:ns,0,0:mpol1)          !z = sum [zmns * sin (mu)], full mesh
      WRITE (51, *) lmnsc(1:ns,0,0:mpol1)          !lam = sum[lmns * sin(mu)], half mesh
!     NOTE: u + lam give a straight magnetic field line
      IF (lasym) THEN
         WRITE (51, *) rmnsc(1:ns,0,0:mpol1)       !r = r+sum [rmns * sin (mu)], full mesh
         WRITE (51, *) zmncc(1:ns,0,0:mpol1)       !z = z+sum [zmnc * cos (mu)], full mesh
         WRITE (51, *) lmncc(1:ns,0,0:mpol1)       !lam = lam+sum[lmnc * cos(mu)], half mesh
      END IF
      WRITE (51, *) chi(1:ns)                      !pol flux, full mesh (included 2*pi factor)
      WRITE (51, *) fpsi(1:ns)                     !R*BT, full mesh
      WRITE (51, *) presf(1:ns)/mu0                !pressure, full mesh (MKS units)
      WRITE (51, *) 1/iotaf(1:ns)                  !q, full mesh

      CLOSE (unit=51)

      DEALLOCATE (rmncc, zmnsc, lmnsc, stat=istat)
      IF (lasym) DEALLOCATE (rmnsc, zmncc, lmncc, stat=istat)

      END SUBROUTINE write_dcon
