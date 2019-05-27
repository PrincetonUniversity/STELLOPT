!-----------------------------------------------------------------------
!     Module:        beams3d_input_mod
!     Authors:       S. Lazerson (lazerson@pppl.gov) M. McMillan (matthew.mcmillan@my.wheaton.edu)
!     Date:          06/20/2012
!     Description:   This module contains the FIELDLINES input namelist and
!                    subroutine which initializes and reads the
!                    FIELDLINES input namelist.
!-----------------------------------------------------------------------
      MODULE beams3d_input_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE beams3d_runtime
      USE beams3d_lines, ONLY: nparticles
      USE beams3d_grid, ONLY: nr, nphi, nz, rmin, rmax, zmin, zmax, &
                              phimin, phimax, vc_adapt_tol, nte, nne, nti,&
                              nzeff, npot
      USE safe_open_mod, ONLY: safe_open
      USE mpi_params

!-----------------------------------------------------------------------
!     Module Variables
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      INCLUDE 'mpif.h'                                                          ! MPI
!DEC$ ENDIF  

!-----------------------------------------------------------------------
!     Input Namelists
!         &beams3d_input
!            nr             Number of radial gridpoints
!            nphi           Number of toroidal gridpoints
!            nz             Number of vertical gridpoints
!            rmin           Minimum radial extent of grid [m]
!            rmax           Maximum radial extent of grid [m]
!            phimin         Minimum toroidal extent of grid [radians]
!            phimax         Maximum toroidal extent of grid [radians]
!            zmin           Minimum vertical extent of grid [m]
!            zmax           Maximum vertical extent of grid [m]
!            mu             Diffusion coefficient
!            r_start        Radial starting locations for fieldlines [m]
!            phi_start      Toroidal starting locations for fieldlines [radians]
!            phi_end        Toroidal ending locations for fieldlines [radians]
!            z_start        Vertical starting locations for fieldlines [m]
!            npoinc         Number of points (per field period) in which data is saved
!            dphi           Fieldlines following stepsize [radians]
!            follow_tol     Tollerance for fieldline following (LSODE and NAG)
!            vc_adapt_tol   Tollerance for adaptive integration using Virtual casing
!                           (note set to negative value to use non-adaptive integration)
!            int_type       Field line integration method
!                           'NAG','LSODE','RKH68'
!
!            NOTE:  Some grid parameters may be overriden (such as
!                   phimin and phimax) to properly represent a given
!                   field period.
!-----------------------------------------------------------------------
      NAMELIST /beams3d_input/ nr, nphi, nz, rmin, rmax, zmin, zmax,&
                                  phimin, phimax, nparticles_start,&
                                  r_start_in, phi_start_in, z_start_in,&
                                  vll_start_in,&
                                  npoinc, follow_tol, t_end_in, mu_start_in,&
                                  charge_in, mass_in, Zatom_in, &
                                  vc_adapt_tol, int_type, Adist_beams,&
                                  Asize_beams, Div_beams, E_beams,&
                                  mass_beams, charge_beams, Zatom_beams, r_beams,&
                                  z_beams, phi_beams, TE_AUX_S, TE_AUX_F,&
                                  NE_AUX_S, NE_AUX_F, TI_AUX_S, TI_AUX_F, &
                                  POT_AUX_S, POT_AUX_F, ZEFF_AUX_S, ZEFF_AUX_F, &
                                  P_beams, ldebug, ne_scale, te_scale, ti_scale, &
                                  zeff_scale
      
!-----------------------------------------------------------------------
!     Subroutines
!         read_beams3d_input:   Reads beams3d_input namelist
!-----------------------------------------------------------------------
      CONTAINS
      
      SUBROUTINE read_beams3d_input(filename, istat)
      CHARACTER(*), INTENT(in) :: filename
      INTEGER, INTENT(out) :: istat
      LOGICAL :: lexist
      INTEGER :: iunit, local_master, i1
      ! Initializations
      local_master = 0
      nr     = 101
      nphi   = 360
      nz     = 101
      rmin   =  0.0_rprec
      rmax   =  1.0_rprec
      zmin   = -1.0_rprec
      zmax   =  1.0_rprec
      phimin =  0.0_rprec
      phimax =  pi2
      nparticles_start = 10

      r_start_in   = -1
      z_start_in   = -1
      phi_start_in = -1
      vll_start_in = -1
      t_end_in     = -1
      mu_start_in  = -1
      mass_in      = -1
      charge_in    = -1
      Zatom_in     = -1

      Adist_beams = 1.0_rprec
      Asize_beams = -1.0_rprec
      Div_beams = 1.0_rprec
      r_beams = 1.0_rprec
      z_beams = 0.0_rprec
      phi_beams = 0.0_rprec
      E_beams = 0.0_rprec
      mass_beams = 1.0_rprec
      charge_beams = 0.0_rprec
      Zatom_beams = 1.0_rprec
      P_beams = 1.0_rprec
      TE_AUX_S = -1
      TE_AUX_F = -1
      NE_AUX_S = -1
      NE_AUX_F = -1
      TI_AUX_S = -1
      TI_AUX_F = -1
      ZEFF_AUX_S = -1
      ZEFF_AUX_F = -1
      POT_AUX_S = -1
      POT_AUX_F = -1
      npoinc = 1
      follow_tol   = 1.0D-7
      vc_adapt_tol = 1.0D-5
      int_type = "LSODE"
      ldebug = .false.
      ne_scale = 1.0
      te_scale = 1.0
      ti_scale = 1.0
      zeff_scale = 1.0
      ! Read namelist
!      IF (ithread == local_master) THEN
         istat=0
         iunit=12
         INQUIRE(FILE=TRIM(filename),EXIST=lexist)
         IF (.not.lexist) stop 'Could not find input file'
         CALL safe_open(iunit,istat,TRIM(filename),'old','formatted')
         IF (istat /= 0) CALL handle_err(NAMELIST_READ_ERR,'beams3d_input in: input.'//TRIM(id_string),istat)
         READ(iunit,NML=beams3d_input,IOSTAT=istat)
         IF (istat /= 0) CALL handle_err(NAMELIST_READ_ERR,'beams3d_input in: input.'//TRIM(id_string),istat)
         CLOSE(iunit)
         NE_AUX_F = NE_AUX_F*ne_scale
         TE_AUX_F = TE_AUX_F*te_scale
         TI_AUX_F = TI_AUX_F*ti_scale
         ZEFF_AUX_F = ZEFF_AUX_F*zeff_scale
         lbeam = .true.
         IF (r_start_in(1) /= -1) lbeam = .false.
         IF (lbeam) lcollision = .true.
         nbeams = 0
         DO WHILE ((Asize_beams(nbeams+1) >= 0.0).and.(nbeams<MAXBEAMS))
            nbeams = nbeams + 1
         END DO
         nte = 0
         DO WHILE ((TE_AUX_S(nte+1) >= 0.0).and.(nte<MAXPROFLEN))
            nte = nte + 1
         END DO
         nne = 0
         DO WHILE ((NE_AUX_S(nne+1) >= 0.0).and.(nne<MAXPROFLEN))
            nne = nne + 1
         END DO
         nti = 0
         DO WHILE ((TI_AUX_S(nti+1) >= 0.0).and.(nti<MAXPROFLEN))
            nti = nti + 1
         END DO
         nzeff = 0
         DO WHILE ((ZEFF_AUX_S(nti+1) >= 0.0).and.(nzeff<MAXPROFLEN))
            nzeff = nzeff + 1
         END DO
         npot = 0
         DO WHILE ((POT_AUX_S(nti+1) >= 0.0).and.(npot<MAXPROFLEN))
            npot = npot + 1
         END DO

         nparticles = 0
         DO WHILE ((r_start_in(nparticles+1) >= 0.0).and.(nparticles<MAXPARTICLES))
            nparticles = nparticles + 1
         END DO
!      END IF

!DEC$ IF DEFINED (HDF5_PAR)
      ! Makes sure that NPARTICLES is divisible by the number of processes
      ! Needed for HDF5 parallel writes.
      IF (lbeam) THEN
         i1 = nparticles_start/nprocs_beams
         IF (i1*nprocs_beams .ne. nparticles_start) THEN
            nparticles_start = (i1+1)*nprocs_beams
         END IF
      END IF
!DEC$ ENDIF

      END SUBROUTINE read_beams3d_input

      SUBROUTINE write_beams3d_namelist(iunit_out, istat)
      INTEGER, INTENT(in) :: iunit_out
      INTEGER, INTENT(out) :: istat
      INTEGER :: ik, n
      CHARACTER(LEN=*), PARAMETER :: outboo  = "(2X,A,1X,'=',1X,L1)"
      CHARACTER(LEN=*), PARAMETER :: outint  = "(2X,A,1X,'=',1X,I0)"
      CHARACTER(LEN=*), PARAMETER :: outflt  = "(2X,A,1X,'=',1X,ES22.12E3)"
      CHARACTER(LEN=*), PARAMETER :: outexp  = "(2X,A,1X,'=',1X,ES22.12E3)"
      CHARACTER(LEN=*), PARAMETER :: outcmp  = "(2x,A,1X,'=','(',i3,',',i3,')')"
      CHARACTER(LEN=*), PARAMETER :: outstr  = "(2X,A,1X,'=',1X,'''',A,'''')"
      CHARACTER(LEN=*), PARAMETER :: onevar  = "(2X,A,1X,'=',1X,L1,2(2X,A,1X,'=',1X,ES22.12E3))"
      CHARACTER(LEN=*), PARAMETER :: vecvar  = "(2X,A,'(',I3.3,')',1X,'=',1X,ES22.12E3)"
      CHARACTER(LEN=*), PARAMETER :: vecvar2  = "(2X,A,'(',I3.3,',',I3.3,')',1X,'=',1X,ES22.12E3)"
      istat = 0
      WRITE(iunit_out,'(A)') '&BEAMS3D_INPUT'
      WRITE(iunit_out,'(A)') '!---------- General Parameters ------------'
      WRITE(iunit_out,outint) 'NR',nr
      WRITE(iunit_out,outint) 'NZ',nz
      WRITE(iunit_out,outint) 'NPHI',nphi
      WRITE(iunit_out,outflt) 'RMIN',rmin
      WRITE(iunit_out,outflt) 'RMAX',rmax
      WRITE(iunit_out,outflt) 'ZMIN',zmin
      WRITE(iunit_out,outflt) 'ZMAX',zmax
      WRITE(iunit_out,outflt) 'PHIMIN',phimin
      WRITE(iunit_out,outflt) 'PHIMAX',phimax
      WRITE(iunit_out,outint) 'NPOINC',npoinc
      WRITE(iunit_out,outstr) 'INT_TYPE',TRIM(int_type)
      WRITE(iunit_out,outflt) 'FOLLOW_TOL',follow_tol
      WRITE(iunit_out,outflt) 'VC_ADAPT_TOL',vc_adapt_tol
      WRITE(iunit_out,outint) 'NPARTICLES_START',nparticles_start
      IF (lbeam) THEN
         WRITE(iunit_out,"(A)") '!---------- Profiles ------------'
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'NE_AUX_S',(ne_aux_s(n), n=1,nne)
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'NE_AUX_F',(ne_aux_f(n), n=1,nne)
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'TE_AUX_S',(te_aux_s(n), n=1,nte)
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'TE_AUX_F',(te_aux_f(n), n=1,nte)
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'TI_AUX_S',(ti_aux_s(n), n=1,nti)
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'TI_AUX_F',(ti_aux_f(n), n=1,nti)
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'ZEFF_AUX_S',(zeff_aux_s(n), n=1,nzeff)
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'ZEFF_AUX_F',(zeff_aux_f(n), n=1,nzeff)
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'POT_AUX_S',(zeff_aux_s(n), n=1,npot)
         WRITE(iunit_out,"(2X,A,1X,'=',4(1X,ES22.12E3))") 'POT_AUX_F',(zeff_aux_f(n), n=1,npot)
         DO n = 1, nbeams
            WRITE(iunit_out,"(A,I2.2)") '!---- BEAM #',n
            WRITE(iunit_out,vecvar) 'T_END_IN',n,t_end_in(n)
            WRITE(iunit_out,vecvar) 'DIV_BEAMS',n,div_beams(n)
            WRITE(iunit_out,vecvar) 'ADIST_BEAMS',n,adist_beams(n)
            WRITE(iunit_out,vecvar) 'ASIZE_BEAMS',n,asize_beams(n)
            WRITE(iunit_out,vecvar) 'MASS_BEAMS',n,mass_beams(n)
            WRITE(iunit_out,vecvar) 'ZATOM_BEAMS',n,zatom_beams(n)
            WRITE(iunit_out,vecvar) 'CHARGE_BEAMS',n,charge_beams(n)
            WRITE(iunit_out,vecvar) 'E_BEAMS',n,e_beams(n)
            WRITE(iunit_out,vecvar) 'P_BEAMS',n,p_beams(n)
            WRITE(iunit_out,vecvar2) 'R_BEAMS',n,1,r_beams(n,1)
            WRITE(iunit_out,vecvar2) 'PHI_BEAMS',n,1,phi_beams(n,1)
            WRITE(iunit_out,vecvar2) 'Z_BEAMS',n,1,z_beams(n,1)
            WRITE(iunit_out,vecvar2) 'R_BEAMS',n,2,r_beams(n,2)
            WRITE(iunit_out,vecvar2) 'PHI_BEAMS',n,2,phi_beams(n,2)
            WRITE(iunit_out,vecvar2) 'Z_BEAMS',n,2,z_beams(n,2)
         END DO
      ELSE
         n = COUNT(t_end_in > -1)
         WRITE(iunit_out,"(2X,A,1X,'=',I0,'*',ES22.12E3)") 'T_END_IN',n,MAXVAL(t_end_in)
      END IF
      WRITE(iunit_out,'(A)') '/'

      END SUBROUTINE write_beams3d_namelist

      SUBROUTINE BCAST_BEAMS3D_INPUT(local_master,comm,istat)
!DEC$ IF DEFINED (MPI_OPT)
      USE mpi
!DEC$ ENDIF
      IMPLICIT NONE
      INTEGER, INTENT(inout) :: comm
      INTEGER, INTENT(in)    :: local_master
      INTEGER, INTENT(inout) :: istat
      IF (istat .ne. 0) RETURN
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BCAST(lbeam, 1, MPI_LOGICAL, local_master, comm,istat)
      CALL MPI_BCAST(nr,1,MPI_INTEGER, local_master, comm,istat)
      CALL MPI_BCAST(nphi,1,MPI_INTEGER, local_master, comm,istat)
      CALL MPI_BCAST(nz,1,MPI_INTEGER, local_master, comm,istat)

      CALL MPI_BCAST(nbeams,1,MPI_INTEGER, local_master, comm,istat)
      CALL MPI_BCAST(nparticles_start,1,MPI_INTEGER, local_master, comm,istat)
      CALL MPI_BCAST(nparticles,1,MPI_INTEGER, local_master, comm,istat)
      CALL MPI_BCAST(npoinc,1,MPI_INTEGER, local_master, comm,istat)
      CALL MPI_BCAST(rmin,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(rmax,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(zmin,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(zmax,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(phimin,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(phimax,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(vc_adapt_tol,1,MPI_REAL8, local_master, comm,istat)

      CALL MPI_BCAST(nte,1,MPI_INTEGER, local_master, comm,istat)
      CALL MPI_BCAST(nne,1,MPI_INTEGER, local_master, comm,istat)
      CALL MPI_BCAST(nti,1,MPI_INTEGER, local_master, comm,istat)
      CALL MPI_BCAST(TE_AUX_S,MAXPROFLEN,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(TE_AUX_F,MAXPROFLEN,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(NE_AUX_S,MAXPROFLEN,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(NE_AUX_F,MAXPROFLEN,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(TI_AUX_S,MAXPROFLEN,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(TI_AUX_F,MAXPROFLEN,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(te_scale,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(ne_scale,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(ti_scale,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(zeff_scale,1,MPI_REAL8, local_master, comm,istat)

      CALL MPI_BCAST(t_end_in,MAXPARTICLES,MPI_REAL8, local_master, comm,istat)

      IF (lbeam) THEN
          CALL MPI_BCAST(Adist_beams,MAXBEAMS,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(Asize_beams,MAXBEAMS,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(Div_beams,MAXBEAMS,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(E_beams,MAXBEAMS,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(r_beams,MAXBEAMS*2,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(z_beams,MAXBEAMS*2,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(phi_beams,MAXBEAMS*2,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(mass_beams,MAXBEAMS,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(charge_beams,MAXBEAMS,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(Zatom_beams,MAXBEAMS,MPI_REAL8, local_master, comm,istat)
      ELSE
          CALL MPI_BCAST(r_start_in,nparticles,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(z_start_in,nparticles,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(phi_start_in,nparticles,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(mu_start_in,nparticles,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(vll_start_in,nparticles,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(mass_in,nparticles,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(charge_in,nparticles,MPI_REAL8, local_master, comm,istat)
          CALL MPI_BCAST(Zatom_in,nparticles,MPI_REAL8, local_master, comm,istat)
      END IF

      CALL MPI_BCAST(follow_tol,1,MPI_REAL8, local_master, comm,istat)
      CALL MPI_BCAST(int_type, 256, MPI_CHARACTER, local_master, comm,istat)
!DEC$ ENDIF
      END SUBROUTINE BCAST_BEAMS3D_INPUT

      END MODULE beams3d_input_mod
