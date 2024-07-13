      SUBROUTINE allocate_saddle_coils
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE boundary, ONLY: nfp
      USE saddle_coils
      USE Vcoilpts
      USE Vwire
      USE coils
      USE mpi_params                                         !mpi stuff
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, n, nc, ncoils, modes, status
!-----------------------------------------------

!     ALLOCATE saddle coil arrays

      ALLOCATE(saddle(nsad_coils_per_period), stat=status)
      IF (myid.EQ.master .AND. status.NE.0) 
     1   STOP "Cannot allocate saddle coils"


!     ALLOCATE saddle coil x,y,z arrays

      ALLOCATE(x_sad(nwdim1,ncdim,nfdim), y_sad(nwdim1,ncdim,nfdim),
     1         z_sad(nwdim1,ncdim,nfdim), stat=status)
      IF (myid.EQ.master .AND. status.NE.0) 
     1   STOP "Cannot allocate saddle xyz arrays"

      ALLOCATE(u_sad(nwdim1,ncdim), v_sad(nwdim1,ncdim),
     1         stat=status)
      IF (myid.EQ.master .AND. status.NE.0) 
     1   STOP "Cannot allocate saddle u-v array"

      nc = nsad_coils_per_period
      IF (myid .EQ. master) THEN 
        IF (nc .le. 0) STOP 'nsad_coils_per_period must > 0!'
      END IF

!     Number of coil types
      nsad_coils = nc * nfp
      IF (myid .EQ. master) THEN
        IF (nsad_coils .gt. ncdim) STOP 'nsad_coils > ncdim'
      END IF

      nsmid = nc/2
      nsodd = MOD(nc,2)
      IF (myid .EQ. master) THEN
        IF (nsodd .EQ. 1) STOP 'nsad_coils_per_period should be even'
      END IF

!     Number of unique coil currents
      num_cursad = 0
      IF (myid .EQ. master) THEN 
        DO i = 1, nsmid
           IF (nsad_group(i) .le. 0) STOP 'nsad_group .le. 0'
        END DO
      END IF

      DO i = 1,nsmid
!     Allocate only those coil components which are needed. There
!     is no need to store the Fourier coefficients for coils whose
!     locations in real space will be computed from stellarator
!     symmetry.

         ALLOCATE(saddle(i)%v_c(0:nsad_v),
     1            saddle(i)%v_s(0:nsad_v),
     2            saddle(i)%u_c(0:nsad_u),
     3            saddle(i)%u_s(0:nsad_u),
     4            stat = status)
      IF (myid.EQ.master .AND. status.NE.0) 
     1   STOP "Cannot allocate saddle components"
      END DO               !! DO i = ... loop over half of coils

      ncoils = nsad_coils                !total no. saddle coils
      nc = ncoils / nfp                  !coils per field period
      IF (ncoils .LE. 0) RETURN

!     initialize the variables to values of unique coil parameters

      DO n = 1, nsmid
        sad_phi0(n) = twopi*sad_v0(n)/nfp
        sad_theta0(n) = twopi*sad_u0(n)
      END DO

      DO i = 1, nsmid
         modes = 0
         saddle(i)%v_c(modes) = sad_v_c(i,modes)
         saddle(i)%v_s(modes) = 0
         DO modes = 1,nsad_v
            saddle(i)%v_c(modes) = sad_v_c(i,modes)
            saddle(i)%v_s(modes) = sad_v_s(i,modes)
         END DO

         modes = 0
         saddle(i)%u_c(modes) = 0
         saddle(i)%u_s(modes) = 0
         IF (nsad_u .gt. 0) THEN
            saddle(i)%u_c(modes) = sad_u_c(i,modes)
            DO modes = 1,nsad_u
               saddle(i)%u_c(modes) = sad_u_c(i,modes)
               saddle(i)%u_s(modes) = sad_u_s(i,modes)
            END DO
         END IF
      END DO

      IF (lsaddle .and. lctrlpt .and. (nsad_u .ne. nsad_v)) THEN
         IF (myid .eq. master) THEN
            PRINT *, 'MUST HAVE NSAD_U = NSAD_V IF LCTRLPT = T'
            STOP
         END IF
      END IF
      IF (lsaddle .and. lspline .and. (nsad_v .lt. 7)) THEN
         IF (myid .eq. master) THEN
            PRINT *, 'MUST HAVE NSAD_V >= 7 IF LSADDLE = T'
            STOP
         END IF
      END IF
      IF (lsaddle .and. lspline .and. (nsad_u .lt. 7)) THEN
         IF (myid .eq. master) THEN
            PRINT *, 'MUST HAVE NSAD_U >= 7 IF LSADDLE = T'
            STOP
         END IF
      END IF

!     Set currents from input for each unique current group

      DO i = 1, nsmid
         saddle(i)%current = cursad(nsad_group(i))*csad_scl(i)
      END DO

!     Impose symmetry on saddle coil-coil penalty weights, EXPonents

      DO i = 1,nsmid
         dsc_wgt(nc+1-i) = dsc_wgt(i)
         dsc_exp(nc+1-i) = dsc_EXP(i)
         dsc_tgt(nc+1-i) = dsc_tgt(i)
         rs_wgt(nc+1-i) = rs_wgt(i)
         rs_exp(nc+1-i) = rs_EXP(i)
         rs_tgt(nc+1-i) = rs_tgt(i)
      END DO

      nsad_unique_coils = nsmid

      END SUBROUTINE allocate_saddle_coils
