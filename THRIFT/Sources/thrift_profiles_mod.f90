!-----------------------------------------------------------------------
!     Module:        thrift_profiles_mod
!     Authors:       L. van Ham, S. Lazerson
!     Date:          11/XX/2022
!     Description:   This module contains the variables and subroutines
!                    for handling profile reading and profile lookup. 
!-----------------------------------------------------------------------
MODULE thrift_profiles_mod
    !-------------------------------------------------------------------
    !     Libraries
    !-------------------------------------------------------------------
    USE stel_kinds, ONLY: rprec
    USE thrift_runtime
    !-------------------------------------------------------------------
    !     Module Variables
    !          lverb         Logical to control screen output
    !-------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: nrho_prof, nt_prof, nion_prof
    INTEGER, DIMENSION(:), POINTER :: Zatom_prof
    REAL(rprec) :: rhomin, rhomax, tmin, tmax, eps1, eps2
    REAL(rprec), DIMENSION(:), POINTER :: raxis_prof, taxis_prof, &
                                          Matom_prof, hr, hri, ht, hti
    REAL(rprec), DIMENSION(:,:,:), POINTER :: NE3D, TE3D, P3D
    REAL(rprec), DIMENSION(:,:,:,:), POINTER :: NI4D, TI4D
    INTEGER :: win_raxis_prof, win_taxis_prof, win_NE3D, win_TE3D, &
               win_NI4D, win_TI4D, win_hr, win_ht, win_hri, win_hti, &
               win_Matom_prof, win_Zatom_prof, win_P3D
!-----------------------------------------------------------------------
!     Input Namelists
!         NONE
!-----------------------------------------------------------------------
      
!-----------------------------------------------------------------------
!     Subroutines
!         read_thrift_profh5: Reads and loads the profiles from HDF5
!         get_prof_ne:        Returns the electron density [m^-3]
!         get_prof_te:        Returns the electron temperature [eV]
!         get_prof_ni:        Returns the ion density [m^-3]
!         get_prof_ti:        Returns the ion temperature [eV]
!         get_prof_p:         Returns the plasma pressure [Pa]
!         get_prof_zeff:      Returns the <Zeff>=sum(n*n*Z)/sum(n*n) [-]
!         get_prof_coulln:    Returns the Coulomb Logarithm [-]
!-----------------------------------------------------------------------
      PUBLIC  :: read_thrift_profh5, get_prof_ne, get_prof_te, &
                 get_prof_ni, get_prof_ti, get_prof_p, free_profiles
      PRIVATE :: setup_grids
      CONTAINS

      SUBROUTINE read_thrift_profh5(filename)
      USE mpi_inc
      USE mpi_params
      USE mpi_sharmem
      USE EZspline
      USE EZspline_obj
#if defined(LHDF5)
      USE ez_hdf5
#endif
      IMPLICIT NONE
      CHARACTER(*), INTENT(in) :: filename
      INTEGER :: i, ier
      INTEGER :: bcs0(2)
      TYPE(EZspline2_r8) :: temp_spl2d
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: temp2d, pres2d
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: temp_ni,temp_ti
      bcs0=(/ 0, 0/)
      ierr_mpi = 0
      IF (lverb) THEN
         WRITE(6,'(A)')  '----- Reading Profile File -----'
         WRITE(6,'(A)')  '   FILE: '//TRIM(filename)
      END IF
      IF (myid_sharmem == master) THEN
         CALL open_hdf5(TRIM(filename),fid,ier,LCREATE=.false.)
         IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,TRIM(filename),ier)
         CALL read_scalar_hdf5(fid,'nrho',ier,INTVAR=nrho_prof)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'nrho_prof',ier)
         CALL read_scalar_hdf5(fid,'nt',ier,INTVAR=nt_prof)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'nt_prof',ier)
         CALL read_scalar_hdf5(fid,'nion',ier,INTVAR=nion_prof)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'nion_prof',ier)
      END IF
      CALL MPI_BARRIER(MPI_COMM_SHARMEM,ierr_mpi)
      ! Broadcast the helpers
      CALL MPI_BCAST(nrho_prof,1,MPI_INTEGER,master,MPI_COMM_SHARMEM,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'read_thrift_profh5: nrho_prof',ierr_mpi)
      CALL MPI_BCAST(nt_prof,1,MPI_INTEGER,master,MPI_COMM_SHARMEM,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'read_thrift_profh5: nt_prof',ierr_mpi)
      CALL MPI_BCAST(nion_prof,1,MPI_INTEGER,master,MPI_COMM_SHARMEM,ierr_mpi)
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'read_thrift_profh5: nion_prof',ierr_mpi)
      ! Allocate the shared memory objects
      CALL mpialloc(raxis_prof, nrho_prof, myid_sharmem, 0, MPI_COMM_SHARMEM, win_raxis_prof)
      CALL mpialloc(taxis_prof, nt_prof,   myid_sharmem, 0, MPI_COMM_SHARMEM, win_taxis_prof)
      CALL mpialloc(Zatom_prof, nion_prof,   myid_sharmem, 0, MPI_COMM_SHARMEM, win_Zatom_prof)
      CALL mpialloc(Matom_prof, nion_prof,   myid_sharmem, 0, MPI_COMM_SHARMEM, win_Matom_prof)
      CALL mpialloc(NE3D, 4, nt_prof, nrho_prof, myid_sharmem, 0, MPI_COMM_SHARMEM, win_NE3D)
      CALL mpialloc(TE3D, 4, nt_prof, nrho_prof, myid_sharmem, 0, MPI_COMM_SHARMEM, win_TE3D)
      CALL mpialloc(P3D,  4, nt_prof, nrho_prof, myid_sharmem, 0, MPI_COMM_SHARMEM, win_P3D)
      CALL mpialloc(NI4D, 4, nt_prof, nrho_prof, nion_prof, myid_sharmem, 0, MPI_COMM_SHARMEM, win_NI4D)
      CALL mpialloc(TI4D, 4, nt_prof, nrho_prof, nion_prof, myid_sharmem, 0, MPI_COMM_SHARMEM, win_TI4D)
      IF (myid_sharmem == master) THEN
         ! Get the axis arrays
         CALL read_var_hdf5(fid,'raxis_prof',nrho_prof,ier,DBLVAR=raxis_prof)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'raxis_prof',ier)
         CALL read_var_hdf5(fid,'taxis_prof',nt_prof,ier,DBLVAR=taxis_prof)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'taxis_prof',ier)
         CALL read_var_hdf5(fid,'Z_prof',nt_prof,ier,INTVAR=Zatom_prof)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'Zatom_prof',ier)
         CALL read_var_hdf5(fid,'mass_prof',nt_prof,ier,DBLVAR=Matom_prof)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'Matom_prof',ier)
         IF (lverb) THEN
            WRITE(6,'(A,F9.5,A,F9.5,A,I4)') '   RHO  = [',minval(raxis_prof),',',maxval(raxis_prof),'];  NRHO: ',nrho_prof
            WRITE(6,'(A,F9.5,A,F9.5,A,I4)') '   T    = [',minval(taxis_prof),',',maxval(taxis_prof),'];    NT: ',nt_prof
         END IF
         !Loop over profiles
         ALLOCATE(temp2d(nt_prof,nrho_prof),pres2d(nt_prof,nrho_prof))
         ! Now create the spline arrays
         ! NE
         CALL read_var_hdf5(fid,'ne_prof',nt_prof,nrho_prof,ier,DBLVAR=temp2d)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'ne_prof',ier)
         IF (lverb) WRITE(6,'(A,F9.5,A,F9.5,A)') '   Ne   = [', &
                        MINVAL(temp2d)*1E-20,',',MAXVAL(temp2d)*1E-20,'] E20 m^-3'
         CALL EZspline_init(temp_spl2d,nt_prof,nrho_prof,bcs0,bcs0,ier)
         temp_spl2d%x1          = taxis_prof
         temp_spl2d%x2          = raxis_prof
         temp_spl2d%isHermite   = 1
         CALL EZspline_setup(temp_spl2d,temp2d,ier,EXACT_DIM=.true.)
         NE3D = temp_spl2d%fspl
         CALL EZspline_free(temp_spl2d,ier)
         pres2d = temp2d
         ! TE
         CALL read_var_hdf5(fid,'te_prof',nt_prof,nrho_prof,ier,DBLVAR=temp2d)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'te_prof',ier)
         IF (lverb) WRITE(6,'(A,F9.5,A,F9.5,A)') '   Te   = [', &
                        MINVAL(temp2d)*1E-3,',',MAXVAL(temp2d)*1E-3,'] keV'
         CALL EZspline_init(temp_spl2d,nt_prof,nrho_prof,bcs0,bcs0,ier)
         temp_spl2d%x1          = taxis_prof
         temp_spl2d%x2          = raxis_prof
         temp_spl2d%isHermite   = 1
         CALL EZspline_setup(temp_spl2d,temp2d,ier,EXACT_DIM=.true.)
         TE3D = temp_spl2d%fspl
         CALL EZspline_free(temp_spl2d,ier)
         pres2d = pres2d*temp2d ! Te*ne
         ! Now work on Ion grid
         DEALLOCATE(temp2d)
         ALLOCATE(temp_ni(nion_prof,nt_prof,nrho_prof))
         ALLOCATE(temp_ti(nion_prof,nt_prof,nrho_prof))
         ! NI
         CALL read_var_hdf5(fid,'ni_prof',nion_prof,nt_prof,nrho_prof,ier,DBLVAR=temp_ni)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'ni_prof',ier)
         DO i = 1, nion_prof
            CALL EZspline_init(temp_spl2d,nt_prof,nrho_prof,bcs0,bcs0,ier)
            IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'init: ni_prof',ier)
            temp_spl2d%x1          = taxis_prof
            temp_spl2d%x2          = raxis_prof
            temp_spl2d%isHermite   = 1
            CALL EZspline_setup(temp_spl2d,temp_ni(i,:,:),ier,EXACT_DIM=.true.)
            IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'setup: ni_prof',ier)
            NI4D(:,:,:,i) = temp_spl2d%fspl
            CALL EZspline_free(temp_spl2d,ier)
         END DO
         ! TI
         CALL read_var_hdf5(fid,'ti_prof',nion_prof,nt_prof,nrho_prof,ier,DBLVAR=temp_ti)
         IF (ier /= 0) CALL handle_err(HDF5_READ_ERR,'ti_prof',ier)
         DO i = 1, nion_prof
            CALL EZspline_init(temp_spl2d,nt_prof,nrho_prof,bcs0,bcs0,ier)
            IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'init: ti_prof',ier)
            temp_spl2d%x1          = taxis_prof
            temp_spl2d%x2          = raxis_prof
            temp_spl2d%isHermite   = 1
            CALL EZspline_setup(temp_spl2d,temp_ti(i,:,:),ier,EXACT_DIM=.true.)
            IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'setup: ti_prof',ier)
            TI4D(:,:,:,i) = temp_spl2d%fspl
            CALL EZspline_free(temp_spl2d,ier)
         END DO
         DO i = 1, nion_prof
            IF (lverb) WRITE(6,'(A,I1,A,F9.5,A,F9.5,A,I3,A,I2)') '   Ni(',i,')= [', &
                        MINVAL(temp_ni(:,:,i))*1E-20,',',MAXVAL(temp_ni(i,:,:))*1E-20,'] E20 m^-3;  M: ',&
                        NINT(Matom_prof(i)/1.66053906660E-27),' amu;  Z: ',Zatom_prof(i)
            IF (lverb) WRITE(6,'(A,I1,A,F9.5,A,F9.5,A)') '   Ti(',i,')= [', &
                        MINVAL(temp_ti(:,:,i))*1E-3,',',MAXVAL(temp_ti(i,:,:))*1E-3,'] keV'
         END DO
         pres2d = (pres2d + SUM(temp_ti*temp_ni,1))*e_charge
         ! DEALLOCATE Helpers
         DEALLOCATE(temp_ti)
         DEALLOCATE(temp_ni)
         ! Close the HDF5 file
         CALL close_hdf5(fid,ier)
         IF (ier /= 0) CALL handle_err(HDF5_CLOSE_ERR,TRIM(filename),ier)
         ! PRESSURE Spline
         IF (lverb) WRITE(6,'(A,F9.5,A,F9.5,A)') '   P    = [', &
                        MINVAL(pres2d)*1E-3,',',MAXVAL(pres2d)*1E-3,'] kPa'
         CALL EZspline_init(temp_spl2d,nt_prof,nrho_prof,bcs0,bcs0,ier)
         temp_spl2d%x1          = taxis_prof
         temp_spl2d%x2          = raxis_prof
         temp_spl2d%isHermite   = 1
         CALL EZspline_setup(temp_spl2d,pres2d,ier,EXACT_DIM=.true.)
         P3D = temp_spl2d%fspl
         CALL EZspline_free(temp_spl2d,ier)
         DEALLOCATE(pres2d)
      END IF
      CALL MPI_BARRIER(MPI_COMM_SHARMEM,ierr_mpi)
      CALL setup_grids
      RETURN
      END SUBROUTINE read_thrift_profh5

      SUBROUTINE setup_grids
      USE mpi_inc
      USE mpi_params
      USE mpi_sharmem
      USE EZspline_type
      IMPLICIT NONE
      INTEGER :: i
      REAL*8, parameter :: small = 1.e-10_ezspline_r8
      CALL mpialloc(hr, nrho_prof-1, myid_sharmem, 0, MPI_COMM_SHARMEM, win_hr)
      CALL mpialloc(hri, nrho_prof-1, myid_sharmem, 0, MPI_COMM_SHARMEM, win_hri)
      CALL mpialloc(ht, nt_prof-1, myid_sharmem, 0, MPI_COMM_SHARMEM, win_hr)
      CALL mpialloc(hti, nt_prof-1, myid_sharmem, 0, MPI_COMM_SHARMEM, win_hri)
      IF (myid_sharmem == master) THEN
        FORALL(i = 1:nrho_prof-1) hr(i) = raxis_prof(i+1) - raxis_prof(i)
        FORALL(i = 1:nt_prof-1)   ht(i) = taxis_prof(i+1) - taxis_prof(i)
        hri = one / hr
        hti = one / ht
      END IF
      rhomin = MINVAL(raxis_prof)
      rhomax = MAXVAL(raxis_prof)
      tmin = MINVAL(taxis_prof)
      tmax = MAXVAL(taxis_prof)
      eps1 = (rhomax-rhomin)*small
      eps2 = (tmax-tmin)*small
      RETURN
      END SUBROUTINE setup_grids

      SUBROUTINE get_prof_f(rho_val,t_val,SPL,val)
      IMPLICIT NONE
      REAL(rprec), INTENT(in) :: rho_val
      REAL(rprec), INTENT(in) :: t_val
      REAL(rprec), INTENT(in) :: SPL(4,nt_prof,nrho_prof)
      REAL(rprec), INTENT(out) :: val
      INTEGER :: i, j
      INTEGER, parameter :: ict(4)=(/1,0,0,0/)
      REAL*8  :: xparam, yparam
      REAL*8 :: fval(1)
      val = 0
      IF ((rho_val >= rhomin-eps1) .and. (rho_val <= rhomax+eps1) .and. &
          (t_val   >= tmin-eps2)   .and. (t_val   <= tmax+eps2)) THEN
         i = MIN(MAX(COUNT(taxis_prof < t_val),1),nt_prof-1)
         j = MIN(MAX(COUNT(raxis_prof < rho_val),1),nrho_prof-1)
         xparam = (t_val - taxis_prof(i)) * hti(i)
         yparam = (rho_val   - raxis_prof(j)) * hri(j)
         CALL R8HERM2FCN(ict,1,1,fval,i,j,xparam,yparam,&
                         ht(i),hti(i),hr(j),hri(j),&
                         SPL(1,1,1),nt_prof,nrho_prof)
         val = fval(1)
      END IF
      RETURN
      END SUBROUTINE get_prof_f

      SUBROUTINE get_prof_dfdrho(rho_val,t_val,SPL,val)
      IMPLICIT NONE
      REAL(rprec), INTENT(in) :: rho_val
      REAL(rprec), INTENT(in) :: t_val
      REAL(rprec), INTENT(in) :: SPL(4,nt_prof,nrho_prof)
      REAL(rprec), INTENT(out) :: val
      INTEGER :: i, j
      INTEGER, parameter :: ict(4)=(/0,0,1,0/)
      REAL*8  :: xparam, yparam
      REAL*8 :: fval(1)
      val = 0
      IF ((rho_val >= rhomin-eps1) .and. (rho_val <= rhomax+eps1) .and. &
          (t_val   >= tmin-eps2)   .and. (t_val   <= tmax+eps2)) THEN
         i = MIN(MAX(COUNT(taxis_prof < t_val),1),nt_prof-1)
         j = MIN(MAX(COUNT(raxis_prof < rho_val),1),nrho_prof-1)
         xparam = (t_val - taxis_prof(i)) * hti(i)
         yparam = (rho_val   - raxis_prof(j)) * hri(j)
         CALL R8HERM2FCN(ict,1,1,fval,i,j,xparam,yparam,&
                         ht(i),hti(i),hr(j),hri(j),&
                         SPL(1,1,1),nt_prof,nrho_prof)
         val = fval(1)
      END IF
      RETURN
      END SUBROUTINE get_prof_dfdrho

      SUBROUTINE get_prof_dfdt(rho_val,t_val,SPL,val)
      IMPLICIT NONE
      REAL(rprec), INTENT(in) :: rho_val
      REAL(rprec), INTENT(in) :: t_val
      REAL(rprec), INTENT(in) :: SPL(4,nt_prof,nrho_prof)
      REAL(rprec), INTENT(out) :: val
      INTEGER :: i, j
      INTEGER, parameter :: ict(4)=(/0,1,0,0/)
      REAL*8  :: xparam, yparam
      REAL*8 :: fval(1)
      val = 0
      IF ((rho_val >= rhomin-eps1) .and. (rho_val <= rhomax+eps1) .and. &
          (t_val   >= tmin-eps2)   .and. (t_val   <= tmax+eps2)) THEN
         i = MIN(MAX(COUNT(taxis_prof < t_val),1),nt_prof-1)
         j = MIN(MAX(COUNT(raxis_prof < rho_val),1),nrho_prof-1)
         xparam = (t_val - taxis_prof(i)) * hti(i)
         yparam = (rho_val   - raxis_prof(j)) * hri(j)
         CALL R8HERM2FCN(ict,1,1,fval,i,j,xparam,yparam,&
                         ht(i),hti(i),hr(j),hri(j),&
                         SPL(1,1,1),nt_prof,nrho_prof)
         val = fval(1)
      END IF
      RETURN
      END SUBROUTINE get_prof_dfdt

      SUBROUTINE get_prof_ne(rho_val,t_val,val)
      IMPLICIT NONE
      REAL(rprec), INTENT(in) :: rho_val
      REAL(rprec), INTENT(in) :: t_val
      REAL(rprec), INTENT(out) :: val
      CALL get_prof_f(rho_val,t_val,NE3D,val)
      RETURN
      END SUBROUTINE get_prof_ne

      SUBROUTINE get_prof_neprime(rho_val,t_val,val)
      IMPLICIT NONE
      REAL(rprec), INTENT(in) :: rho_val
      REAL(rprec), INTENT(in) :: t_val
      REAL(rprec), INTENT(out) :: val
      CALL get_prof_dfdrho(rho_val,t_val,NE3D,val)
      RETURN
      END SUBROUTINE get_prof_neprime

      SUBROUTINE get_prof_te(rho_val,t_val,val)
      IMPLICIT NONE
      REAL(rprec), INTENT(in) :: rho_val
      REAL(rprec), INTENT(in) :: t_val
      REAL(rprec), INTENT(out) :: val
      CALL get_prof_f(rho_val,t_val,TE3D,val)
      RETURN
      END SUBROUTINE get_prof_te

      SUBROUTINE get_prof_teprime(rho_val,t_val,val)
      IMPLICIT NONE
      REAL(rprec), INTENT(in) :: rho_val
      REAL(rprec), INTENT(in) :: t_val
      REAL(rprec), INTENT(out) :: val
      CALL get_prof_dfdrho(rho_val,t_val,TE3D,val)
      RETURN
      END SUBROUTINE get_prof_teprime

      SUBROUTINE get_prof_ni(rho_val,t_val,iion,val)
      IMPLICIT NONE
      REAL(rprec), INTENT(in) :: rho_val
      REAL(rprec), INTENT(in) :: t_val
      INTEGER,     INTENT(in) :: iion
      REAL(rprec), INTENT(out) :: val
      val = 0
      IF (iion <= nion_prof) CALL get_prof_f(rho_val,t_val,NI4D(:,:,:,iion),val)
      RETURN
      END SUBROUTINE get_prof_ni

      SUBROUTINE get_prof_niprime(rho_val,t_val,iion,val)
      IMPLICIT NONE
      REAL(rprec), INTENT(in) :: rho_val
      REAL(rprec), INTENT(in) :: t_val
      INTEGER,     INTENT(in) :: iion
      REAL(rprec), INTENT(out) :: val
      val = 0
      IF (iion <= nion_prof) CALL get_prof_dfdrho(rho_val,t_val,NI4D(:,:,:,iion),val)
      RETURN
      END SUBROUTINE get_prof_niprime

      SUBROUTINE get_prof_ti(rho_val,t_val,iion,val)
      IMPLICIT NONE
      REAL(rprec), INTENT(in) :: rho_val
      REAL(rprec), INTENT(in) :: t_val
      INTEGER,     INTENT(in) :: iion
      REAL(rprec), INTENT(out) :: val
      val = 0
      IF (iion <= nion_prof) CALL get_prof_f(rho_val,t_val,TI4D(:,:,:,iion),val)
      RETURN
      END SUBROUTINE get_prof_ti

      SUBROUTINE get_prof_tiprime(rho_val,t_val,iion,val)
      IMPLICIT NONE
      REAL(rprec), INTENT(in) :: rho_val
      REAL(rprec), INTENT(in) :: t_val
      INTEGER,     INTENT(in) :: iion
      REAL(rprec), INTENT(out) :: val
      val = 0
      IF (iion <= nion_prof) CALL get_prof_dfdrho(rho_val,t_val,TI4D(:,:,:,iion),val)
      RETURN
      END SUBROUTINE get_prof_tiprime

      SUBROUTINE get_prof_p(rho_val,t_val,val)
      IMPLICIT NONE
      REAL(rprec), INTENT(in) :: rho_val
      REAL(rprec), INTENT(in) :: t_val
      REAL(rprec), INTENT(out) :: val
      INTEGER     :: i
      REAL(rprec) :: nk, tk
      val = 0
      CALL get_prof_f(rho_val,t_val,P3D,val)
      RETURN
      END SUBROUTINE get_prof_p

      SUBROUTINE get_prof_pprime(rho_val,t_val,val)
      IMPLICIT NONE
      REAL(rprec), INTENT(in) :: rho_val
      REAL(rprec), INTENT(in) :: t_val
      REAL(rprec), INTENT(out) :: val
      INTEGER     :: i
      REAL(rprec) :: nk, tk, dn, dt
      val = 0
      CALL get_prof_dfdrho(rho_val,t_val,P3D,val)
      RETURN
      END SUBROUTINE get_prof_pprime

      SUBROUTINE get_prof_zeff(rho_val,t_val,val)
      IMPLICIT NONE
      REAL(rprec), INTENT(in) :: rho_val
      REAL(rprec), INTENT(in) :: t_val
      REAL(rprec), INTENT(out) :: val
      INTEGER     :: i
      REAL(rprec) :: f_top, f_bot, nk
      val = 0; f_top = 0; f_bot = 0
      DO i = 1, nion_prof
         CALL get_prof_ni(rho_val,t_val,i,nk)
         f_top = f_top + nk*nk*Zatom_prof(i)
         f_bot = f_bot + nk*nk
      END DO
      IF (f_bot > 0) val = f_top/f_bot
      RETURN
      END SUBROUTINE get_prof_zeff

      SUBROUTINE get_prof_coulln(rho_val,t_val,val)
      IMPLICIT NONE
      REAL(rprec), INTENT(in) :: rho_val
      REAL(rprec), INTENT(in) :: t_val
      REAL(rprec), INTENT(out) :: val
      INTEGER     :: i
      REAL(rprec) :: ne, te
      val = 0
      ! NRL19 pg34 electron-ion
      CALL get_prof_ne(rho_val,t_val,ne)
      CALL get_prof_te(rho_val,t_val,te)
      val = 24 - log( sqrt(ne*1E-6)/(te) )
      RETURN
      END SUBROUTINE get_prof_coulln

      SUBROUTINE get_prof_etaperp(rho_val,t_val,val)
      IMPLICIT NONE
      REAL(rprec), INTENT(in) :: rho_val
      REAL(rprec), INTENT(in) :: t_val
      REAL(rprec), INTENT(out) :: val
      INTEGER     :: i
      REAL(rprec) :: clog, te, zeff
      val = 0
      ! https://en.wikipedia.org/wiki/Spitzer_resistivity
      CALL get_prof_te(rho_val,t_val,te)
      CALL get_prof_zeff(rho_val,t_val,zeff)
      CALL get_prof_coulln(rho_val,t_val,clog)
      val = 1.0313621201E-04*zeff*clog/(te**1.5)
      RETURN
      END SUBROUTINE get_prof_etaperp

      SUBROUTINE get_prof_etapara(rho_val,t_val,val)
      IMPLICIT NONE
      REAL(rprec), INTENT(in) :: rho_val
      REAL(rprec), INTENT(in) :: t_val
      REAL(rprec), INTENT(out) :: val
      INTEGER     :: i
      REAL(rprec) :: zeff, f
      val = 0
      ! https://en.wikipedia.org/wiki/Spitzer_resistivity
      CALL get_prof_zeff(rho_val,t_val,zeff)
      CALL get_prof_etaperp(rho_val,t_val,val)
      F = (1+1.198*zeff+0.222*zeff*zeff)/(1+2.966*zeff+0.752*zeff*zeff)
      val = val * F
      RETURN
      END SUBROUTINE get_prof_etapara

      SUBROUTINE free_profiles
      USE mpi_sharmem
      IMPLICIT NONE
      IF (ASSOCIATED(raxis_prof)) CALL mpidealloc(raxis_prof,win_raxis_prof)
      IF (ASSOCIATED(taxis_prof)) CALL mpidealloc(taxis_prof,win_taxis_prof)
      IF (ASSOCIATED(Matom_prof)) CALL mpidealloc(Matom_prof,win_Matom_prof)
      IF (ASSOCIATED(Zatom_prof)) CALL mpidealloc(Zatom_prof,win_Zatom_prof)
      IF (ASSOCIATED(hr))         CALL mpidealloc(hr,win_hr)
      IF (ASSOCIATED(hri))        CALL mpidealloc(hri,win_hri)
      IF (ASSOCIATED(ht))         CALL mpidealloc(ht,win_ht)
      IF (ASSOCIATED(hti))        CALL mpidealloc(hti,win_hti)
      IF (ASSOCIATED(NE3D))       CALL mpidealloc(NE3D,win_NE3D)
      IF (ASSOCIATED(TE3D))       CALL mpidealloc(TE3D,win_TE3D)
      IF (ASSOCIATED(P3D))        CALL mpidealloc(TE3D,win_P3D)
      IF (ASSOCIATED(NI4D))       CALL mpidealloc(NI4D,win_NI4D)
      IF (ASSOCIATED(TI4D))       CALL mpidealloc(TI4D,win_TI4D)
      END SUBROUTINE free_profiles

END MODULE thrift_profiles_mod