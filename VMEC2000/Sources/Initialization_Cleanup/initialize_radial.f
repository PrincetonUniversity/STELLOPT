      SUBROUTINE initialize_radial(nsval, ns_old, delt0,
     1                             lscreen, reset_file_name)
      USE vmec_main
      USE vmec_params, ONLY: ntmax 
      USE realspace
      USE vsvd
      USE xstuff
#ifdef _HBANGLE
      USE angle_constraints, ONLY: getrz, store_init_array
#endif
      USE parallel_include_module
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(IN) :: nsval
      INTEGER, INTENT(INOUT) :: ns_old
      CHARACTER(LEN=*), OPTIONAL :: reset_file_name
      REAL(dp), INTENT(OUT) :: delt0
      LOGICAL, INTENT(IN) :: lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: neqs2_old=0
      LOGICAL :: lreset_internal, linterp
      INTEGER :: nsmin, nsmax, i, j, k, l, lk
      REAL(dp) :: tiniton, tinitoff
C-----------------------------------------------
!
!     Allocates memory for radial arrays and initializes radial profiles
!     Loads data (if available) from a reset file
!
C-----------------------------------------------
!
!                INDEX OF LOCAL VARIABLES
!
!        hs      radial mesh size increment
!        irzloff offset in xc array between R,Z,L components
!        neqs    total number of equations to evolve (size of xc)
C-----------------------------------------------
#if defined(SKS)
      CALL second0(tiniton)
#endif
!     Set timestep control parameters
      fsq     = one
      iter2 = 1
      iter1 = iter2
      ijacob = 0
      irst = 1
      res0 = -1
!
!       INITIALIZE MESH-DEPENDENT SCALARS
!
      ns = nsval
      ns1 = ns-1
      delt0 = delt
      hs = one/ns1
      ohs = one/hs
      mns = ns*mnsize
      irzloff = ntmax*mns
      nrzt = nznt*ns
      neqs = 3*irzloff
      neqs1 = neqs + 1
      neqs2 = neqs1 + 1


      IF (grank .EQ. 0) THEN
        WRITE (nthreed, 10) ns, mnmax, ftolv, niter
        IF (lscreen) PRINT 10, ns, mnmax, ftolv, niter
   10 FORMAT(/'  NS = ',i4,' NO. FOURIER MODES = ',i4,' FTOLV = ',
     1   1p,e10.3,' NITER = ',i6)
        IF (lactive) THEN
           IF (lfreeb) THEN
              WRITE(nthreed, 14) nranks, vnranks
              IF (lscreen) PRINT 14, nranks, vnranks
           ELSE
              WRITE (nthreed, 12) nranks
              IF (lscreen) PRINT 12, nranks
           END IF
  12  FORMAT ('  PROCESSOR COUNT - RADIAL: ',i4)
  14  FORMAT ('  PROCESSOR COUNT - RADIAL: ',i4,'  VACUUM: ',i4)
        END IF

      IF (imovephi .gt. 0 .and. lscreen) PRINT *,
     1   'Turning on Ip Matching by varying Phi-Edge.'
      END IF

!
!     ALLOCATE NS-DEPENDENT ARRAYS
!
      lreset_internal = .true.
      linterp = (ns_old.LT.ns .AND. ns_old.NE.0)
      IF (ns_old .EQ. ns) RETURN
      CALL allocate_ns(linterp, neqs2_old)
!
!     SAVE THIS FOR INTERPOLATION
!
      IF (neqs2_old.gt.0 .and. linterp) THEN
#ifdef _HBANGLE
         ns = ns_old
         CALL getrz(xstore)
         ns = ns1+1
#endif
         IF(PARVMEC) THEN
#if defined(MPI_OPT)
           pgc(1:neqs2_old)=pscalxc(1:neqs2_old)*pxstore(1:neqs2_old)
           IF (lfreeb) 
     1        CALL MPI_Bcast(rbsq,SIZE(rbsq),MPI_REAL8,
     2                       0,NS_COMM,MPI_ERR)
#endif
         ELSE
           gc(1:neqs2_old)=scalxc(1:neqs2_old)*xstore(1:neqs2_old)
         END IF
      END IF

!
!     COMPUTE INITIAL R, Z AND MAGNETIC FLUX PROFILES
!

      IF (PARVMEC) THEN
        CALL profil1d_par (pxc, pxcdot, lreset_internal)
      ELSE
        CALL profil1d (xc, xcdot, lreset_internal)
      END IF
      IF (PRESENT(reset_file_name)) THEN
         IF (LEN_TRIM(reset_file_name) .ne. 0) THEN
            CALL load_xc_from_wout(xc(1), xc(1+irzloff),
     &                             xc(1+2*irzloff), lreset_internal,
     &                             ntor, mpol1, ns, reset_file_name)
            IF (PARVMEC) THEN
               CALL Serial2Parallel4X(xc,pxc)
            END IF
         END IF
      END IF

      IF (PARVMEC) THEN
        CALL profil3d_par(pxc(1),pxc(1+irzloff),lreset_internal,linterp)
      ELSE
        CALL profil3d (xc(1), xc(1+irzloff), lreset_internal, linterp)
      END IF
!
!     INTERPOLATE FROM COARSE (ns_old) TO NEXT FINER (ns) RADIAL GRID
!
      IF (linterp) THEN
        IF(PARVMEC) THEN 
          CALL interp_par (pxc, pgc, pscalxc, ns, ns_old)
        ELSE
          CALL interp (xc, gc, scalxc, ns, ns_old)
        END IF
#ifdef _HBANGLE
        CALL store_init_array(xc)
#endif
      END IF

!SPH 012417: move this AFTER interpolation call
      irst = 1
      CALL restart_iter(delt)

      ns_old = ns
      neqs2_old = neqs2

#if defined (SKS)
      CALL second0(tinitoff)
      init_radial_time = init_radial_time + (tinitoff-tiniton)
#endif

      END SUBROUTINE initialize_radial
