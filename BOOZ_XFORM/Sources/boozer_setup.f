      SUBROUTINE boozer_setup
      USE booz_params
      USE booz_persistent
      USE mpi_params, ONLY: myid, master
      IMPLICIT NONE
C-----------------------------------------------
!
!     One-time initialization of the Boozer angle grids and
!     transformation scale factors.  Previously this work was
!     performed inside boozer_coords during its first invocation
!     (jsurf == 0).  It has been factored out so that it can be
!     executed once on every MPI rank prior to the parallel flux
!     surface loop.
!
      CALL setup_booz (ntorsum, ns, mnmax, ohs, xmb, xnb,
     1   sfull, scl, mboz, nboz, mnboz, nu2_b, nu_boz,
     2   nv_boz, nfp, lasym_b)

      IF (lasym_b) THEN
         nu3_b = nu_boz
      ELSE
         nu3_b = nu2_b                 !!ONLY need top half of theta mesh for symmetric plasma
      END IF

      nunv = nu3_b*nv_boz

      CALL foranl (nu3_b, nv_boz, nfp, nunv, lasym_b)

      IF (lscreen .AND. myid .eq. master) THEN
      	WRITE(6, 50) mboz-1, -nboz, nboz, nu_boz, nv_boz
      END IF
      
  50  FORMAT('  0 <= mboz <= ',i4,3x,i4,' <= nboz <= ',i4,/,
     1       '  nu_boz = ',i5,' nv_boz = ',i5,//,
     2       13x,'OUTBOARD (u=0)',14x,'JS',10x,'INBOARD (u=pi)'
     3       /,77('-')/,'  v     |B|vmec    |B|booz    Error',13x,
     4       '|B|vmec    |B|booz    Error'/)

      END SUBROUTINE boozer_setup
