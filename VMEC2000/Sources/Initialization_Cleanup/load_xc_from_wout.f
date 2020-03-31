      SUBROUTINE load_xc_from_wout(rmn, zmn, lmn, lreset, 
     1    ntor_in, mpol1_in, ns_in, reset_file)
      USE read_wout_mod, ONLY: rmnc, zmns, lmns, rmns, zmnc, lmnc,
     1    xm, xn, ntor, ns,
     2    nfp, mnmax, read_wout_file, read_wout_deallocate
      USE vmec_params, ONLY: mscale, nscale, ntmax,
     1                       rcc, rss, rsc, rcs, zsc, zcs, zcc, zss
      USE vmec_dim, ONLY: mpol1
      USE vparams, ONLY: one, zero, rprec
      USE vmec_input, ONLY: lasym
      USE vmec_main, ONLY: lthreed, p5 => cp5, sp, sm, phipf
      USE parallel_include_module, ONLY: rank
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: ns_in, mpol1_in, ntor_in
      REAL(rprec), DIMENSION(ns_in,0:ntor_in,0:mpol1_in,ntmax),
     1   INTENT(out) :: rmn, zmn, lmn
      LOGICAL, INTENT(out) :: lreset
      CHARACTER(LEN=*) :: reset_file
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ierr, mn, m, n, n1, js
      REAL(rprec) :: t1, t2
      REAL(rprec), ALLOCATABLE :: temp(:,:)
C-----------------------------------------------

!
!     THIS ALLOWS SEQUENTIAL RUNNING OF VMEC FROM THE COMMAND LINE
!     i.e., WHEN VMEC INTERNAL ARRAYS ARE NOT KEPT IN MEMORY (compared to sequence file input)
!     THIS IS THE CASE WHEN VMEC IS CALLED FROM, SAY, THE OPTIMIZATION CODE
!
!     SPH 12-13-11: allow for paths in wout file name (as per Ed Lazarus request)
      CALL read_wout_file (reset_file, ierr)
!      CALL read_wout_file (reset_file(5:), ierr)
      reset_file = " "               !nullify so this routine will not be recalled with present reset_file

      IF (ierr .ne. 0.AND.rank.EQ.0) THEN
         PRINT *,' Error opening/reading wout file in VMEC load_xc!'
         RETURN
      END IF

      IF (ns_in .ne. ns.AND.rank.EQ.0) THEN
         PRINT *, 'ns_in (passed to load_xc) != ns (from reading wout)'
         RETURN
      END IF

      IF (ntor_in  .ne. ntor ) STOP 'ntor_in != ntor in load_xc'
      IF (mpol1_in .ne. mpol1) STOP 'mpol1_in != mpol1 in load_xc'
      IF (nfp .eq. 0) STOP 'nfp = 0 in load_xc'

      lreset = .false.               !Signals profil3d NOT to overwrite axis values

      rmn = zero
      zmn = zero
      lmn = zero

      DO mn = 1, mnmax
         m = NINT(xm(mn))
         n = NINT(xn(mn))/nfp
         n1 = ABS(n)
         t1 = one/(mscale(m)*nscale(n1))
         t2 = t1
         IF (n .lt. 0) t2 = -t2
         IF (n .eq. 0) t2 = zero
         rmn(:ns,n1,m,rcc) = rmn(:ns,n1,m,rcc) + t1*rmnc(mn,:ns)
         zmn(:ns,n1,m,zsc) = zmn(:ns,n1,m,zsc) + t1*zmns(mn,:ns)
         lmn(:ns,n1,m,zsc) = lmn(:ns,n1,m,zsc) + t1*lmns(mn,:ns)
         IF (lthreed) THEN
         rmn(:ns,n1,m,rss) = rmn(:ns,n1,m,rss) + t2*rmnc(mn,:ns)
         zmn(:ns,n1,m,zcs) = zmn(:ns,n1,m,zcs) - t2*zmns(mn,:ns)
         lmn(:ns,n1,m,zcs) = lmn(:ns,n1,m,zcs) - t2*lmns(mn,:ns)
         END IF
         IF (lasym) THEN
         rmn(:ns,n1,m,rsc) = rmn(:ns,n1,m,rsc) + t1*rmns(mn,:ns)
         zmn(:ns,n1,m,zcc) = zmn(:ns,n1,m,zcc) + t1*zmnc(mn,:ns)
         lmn(:ns,n1,m,zcc) = lmn(:ns,n1,m,zcc) + t1*lmnc(mn,:ns)
         IF (lthreed) THEN
         rmn(:ns,n1,m,rcs) = rmn(:ns,n1,m,rcs) - t2*rmns(mn,:ns)
         zmn(:ns,n1,m,zss) = zmn(:ns,n1,m,zss) + t2*zmnc(mn,:ns)
         lmn(:ns,n1,m,zss) = lmn(:ns,n1,m,zss) + t2*lmnc(mn,:ns)
         END IF
         END IF
         IF (m .eq. 0) THEN
            zmn(:ns,n1,m,zsc) = zero
            lmn(:ns,n1,m,zsc) = zero
            IF (lthreed) rmn(:ns,n1,m,rss) = zero
            IF (lasym) THEN
            rmn(:ns,n1,m,rsc) = zero
            IF (lthreed) THEN
            zmn(:ns,n1,m,zss) = zero
            lmn(:ns,n1,m,zss) = zero
            END IF
            END IF
         END IF
      END DO

!
!     CONVERT TO INTERNAL FORM FOR (CONSTRAINED) m=1 MODES
!

      IF (lthreed .or. lasym) ALLOCATE (temp(ns_in,0:ntor_in))
      IF (lthreed) THEN
         temp = rmn(:,:,1,rss)
         rmn(:,:,1,rss) = p5*(temp + zmn(:,:,1,zcs))
         zmn(:,:,1,zcs) = p5*(temp - zmn(:,:,1,zcs))
      END IF
      IF (lasym) THEN
         temp = rmn(:,:,1,rsc)
         rmn(:,:,1,rsc) = p5*(temp + zmn(:,:,1,zcc))
         zmn(:,:,1,zcc) = p5*(temp - zmn(:,:,1,zcc))
      END IF

      IF (ALLOCATED(temp)) DEALLOCATE (temp)

!
!     CONVERT lambda TO INTERNAL FULL MESH REPRESENTATION
!
!     START ITERATION AT JS=1
!
      lmn(1,:,0,:) = lmn(2,:,0,:)
      lmn(1,:,1,:) = 2*lmn(2,:,1,:)/(sm(2) + sp(1))
      lmn(1,:,2:,:) = 0
      
      DO m = 0, mpol1, 2
         DO js = 2, ns
            lmn(js,:,m,:) = 2*lmn(js,:,m,:) - lmn(js-1,:,m,:)
         END DO
      END DO

      DO m = 1, mpol1, 2
         DO js = 2, ns
            lmn(js,:,m,:) = (2*lmn(js,:,m,:) 
     1                    - sp(js-1)*lmn(js-1,:,m,:))/sm(js)
         END DO
      END DO

      DO js = 2, ns
         lmn(js,:,:,:) = phipf(js)*lmn(js,:,:,:)
      END DO

      CALL read_wout_deallocate

      END SUBROUTINE load_xc_from_wout
