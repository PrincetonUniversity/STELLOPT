 !-----------------------------------------------------------------------
!     Module:        beams3d_neutdens
!     Authors:       L. van Ham (lucas.van.ham@ipp.mpg.de)
!     Date:          02/12/2025
!     Description:   This subroutine reads the neutral density on a grid
!                    for the neutralizer and sets up splines.
!-----------------------------------------------------------------------

MODULE beams3d_neutdens
    !-----------------------------------------------------------------------
    !     Libraries
    !-----------------------------------------------------------------------
    USE stel_kinds, ONLY: rprec

#if defined(LHDF5)
    USE ez_hdf5
#endif
    USE safe_open_mod, ONLY: safe_open
    USE mpi_inc
    USE mpi_params
    USE mpi_sharmem
    USE EZspline_obj
    USE EZspline
    USE beams3d_runtime

    INTEGER :: n_u, n_v, n_w
    DOUBLE PRECISION :: neut_x0(3), neut_delta_u, neut_delta_v, neut_delta_w
    DOUBLE PRECISION :: neut_dir_u(3), neut_dir_v(3), neut_dir_w(3)
    INTEGER :: win_neutdens_spl
    REAL(rprec), POINTER :: neutdens_fspl(:,:,:,:)
    REAL(rprec) :: hr_u, hri_u, hr_v, hri_v, hr_w, hri_w
    REAL(rprec), ALLOCATABLE ::neut_grid_u(:),neut_grid_v(:), neut_grid_w(:)
    

CONTAINS

SUBROUTINE beams3d_read_neutdens(filename)
    !-----------------------------------------------------------------------
    !     Subroutine: beams3d_reud_neutdens
    !           Reads a .h5 neutral density file and sets up its spline
    !-----------------------------------------------------------------------
    !     Input Variables
    !          filename     name of file to read (e.g. NI21_density.h5)
    !-----------------------------------------------------------------------
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(in)           :: filename
    !-----------------------------------------------------------------------
    !     Local Variables
    !          ier               Error Flag
    !          iunit             File ID
    !          n_u,v,w           Number of gridpoints in each direction 
    !          neut_grid_u,v,w   Grid of coordinates
    !          neut_density      Density at these coordinats 
    !-----------------------------------------------------------------------
        INTEGER :: ier, iunit
        REAL(rprec), ALLOCATABLE :: neut_density(:,:,:)
        INTEGER :: i, MPI_COMM_MASTERS, ierr_mpi
        INTEGER :: bcs1(2), bcs2(2), bcs3(2)
        TYPE(EZspline3_r8) :: spline_local

    !-----------------------------------------------------------------------
    !     Begin Subroutine
    !-----------------------------------------------------------------------
    IF (lverb) THEN
        WRITE(6,'(A)')  '----- READING NEUTRAL DENSITY FROM FILE -----'
    END IF
#if defined(LHDF5)
    IF (myworkid == master) THEN
        IF (lverb) WRITE(6,'(A)')  '   FILE: '//TRIM(filename)
        CALL open_hdf5(TRIM(filename),fid,ier,LCREATE=.false.)
        IF (ier /= 0) CALL handle_err(HDF5_OPEN_ERR,TRIM(filename),ier)
        
        ! Read dimensions
        CALL read_scalar_hdf5(fid,'/dim/u',ier,INTVAR=n_u)
        CALL read_scalar_hdf5(fid,'/dim/v',ier,INTVAR=n_v)
        CALL read_scalar_hdf5(fid,'/dim/w',ier,INTVAR=n_w)

        ! Allocate
        ALLOCATE(neut_density(n_u,n_v,n_w),neut_grid_u(n_u),neut_grid_v(n_v),neut_grid_w(n_w))

        ! Read data
        CALL read_var_hdf5(fid,'/dir/u',3,ier,DBLVAR=neut_dir_u)
        CALL read_var_hdf5(fid,'/dir/v',3,ier,DBLVAR=neut_dir_v)
        CALL read_var_hdf5(fid,'/dir/w',3,ier,DBLVAR=neut_dir_w)
        CALL read_scalar_hdf5(fid,'/delta/u',ier,DBLVAR=neut_delta_u)
        CALL read_scalar_hdf5(fid,'/delta/v',ier,DBLVAR=neut_delta_v)
        CALL read_scalar_hdf5(fid,'/delta/w',ier,DBLVAR=neut_delta_w)
        CALL read_var_hdf5(fid,'/x0',3,ier,DBLVAR=neut_x0)
        CALL read_var_hdf5(fid,'/data',n_u,n_v,n_w,ier,DBLVAR=neut_density)
        CALL read_var_hdf5(fid,'/grid/u',n_u,ier,DBLVAR=neut_grid_u)
        CALL read_var_hdf5(fid,'/grid/v',n_v,ier,DBLVAR=neut_grid_v)
        CALL read_var_hdf5(fid,'/grid/w',n_w,ier,DBLVAR=neut_grid_w)
        IF (lverb) WRITE(6,'(A,I3,A,I3,A,I3)')     '   Dimensions:  ', n_u, ', ', n_v, ', ', n_w
        IF (lverb) WRITE(6,'(A,E10.3,A,E10.3,A)') '   n [m^-3]  : [',MINVAL(neut_density),',',MAXVAL(neut_density),']'

        ! Close the file
        CALL close_hdf5(fid,ier)

        ! Grid spacings
        hr_u = neut_grid_u(2) - neut_grid_u(1)
        hri_u = 1.0_rprec/hr_u
        hr_v = neut_grid_v(2) - neut_grid_v(1)
        hri_v = 1.0_rprec/hr_v
        hr_w = neut_grid_w(2) - neut_grid_w(1)
        hri_w = 1.0_rprec/hr_w

    END IF

    CALL MPI_BCAST(hr_u, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
    CALL MPI_BCAST(hri_u, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
    CALL MPI_BCAST(hr_v, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
    CALL MPI_BCAST(hri_v, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
    CALL MPI_BCAST(hr_w, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
    CALL MPI_BCAST(hri_w, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
    
    ! Find masters
    i = MPI_UNDEFINED
    IF (myid_sharmem == master) i = 0
    CALL MPI_COMM_SPLIT(MPI_COMM_BEAMS, i, myworkid, MPI_COMM_MASTERS, ierr_mpi)

    ! Broadcast from rank 0 to all other masters
    IF (myid_sharmem == master) THEN
        CALL MPI_BCAST(n_u, 1, MPI_INTEGER, 0, MPI_COMM_MASTERS, ierr_mpi)
        CALL MPI_BCAST(n_v, 1, MPI_INTEGER, 0, MPI_COMM_MASTERS, ierr_mpi)
        CALL MPI_BCAST(n_w, 1, MPI_INTEGER, 0, MPI_COMM_MASTERS, ierr_mpi)
        CALL MPI_BCAST(neut_x0, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_MASTERS, ierr_mpi)
        CALL MPI_BCAST(neut_delta_u, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_MASTERS, ierr_mpi)
        CALL MPI_BCAST(neut_delta_v, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_MASTERS, ierr_mpi)
        CALL MPI_BCAST(neut_delta_w, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_MASTERS, ierr_mpi)
        CALL MPI_BCAST(neut_dir_u, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_MASTERS, ierr_mpi)
        CALL MPI_BCAST(neut_dir_v, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_MASTERS, ierr_mpi)
        CALL MPI_BCAST(neut_dir_w, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_MASTERS, ierr_mpi)

        ! Allocate remaining masters
        if (myworkid /= master) ALLOCATE(neut_density(n_u,n_v,n_w),neut_grid_u(n_u),neut_grid_v(n_v),neut_grid_w(n_w))

        CALL MPI_BCAST(neut_density, n_u*n_v*n_w, MPI_DOUBLE_PRECISION, 0, MPI_COMM_MASTERS, ierr_mpi)
        CALL MPI_BCAST(neut_grid_u, n_u, MPI_DOUBLE_PRECISION, 0, MPI_COMM_MASTERS, ierr_mpi)
        CALL MPI_BCAST(neut_grid_v, n_v, MPI_DOUBLE_PRECISION, 0, MPI_COMM_MASTERS, ierr_mpi)
        CALL MPI_BCAST(neut_grid_w, n_w, MPI_DOUBLE_PRECISION, 0, MPI_COMM_MASTERS, ierr_mpi)
        
        ! Set-up spline on the shared memory window
        ! spline is in u,v,w, all going from [0 1]
        bcs1 = (/ 0, 0 /)
        bcs2 = (/ 0, 0 /)
        bcs3 = (/ 0, 0 /)
        CALL EZspline_init(spline_local,n_u,n_v,n_w,bcs1,bcs2,bcs3,ier)
        IF (ier /= 0) CALL handle_err(EZSPLINE_ERR, 'neutdens_spl_init', ier)
        spline_local%isHermite   = 1
        spline_local%x1   = neut_grid_u
        spline_local%x2   = neut_grid_v
        spline_local%x3   = neut_grid_w
        CALL EZspline_setup(spline_local, neut_density, ier, EXACT_DIM=.TRUE.)
        IF (ier /= 0) CALL handle_err(EZSPLINE_ERR, 'neutdens_spl_setup', ier)
    END IF
    
    ! Broadcast to MPI threads now
    CALL MPI_BCAST(n_u, 1, MPI_INTEGER, 0, MPI_COMM_SHARMEM, ierr_mpi)
    CALL MPI_BCAST(n_v, 1, MPI_INTEGER, 0, MPI_COMM_SHARMEM, ierr_mpi)
    CALL MPI_BCAST(n_w, 1, MPI_INTEGER, 0, MPI_COMM_SHARMEM, ierr_mpi)
    CALL MPI_BCAST(neut_x0, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SHARMEM, ierr_mpi)
    CALL MPI_BCAST(neut_delta_u, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SHARMEM, ierr_mpi)
    CALL MPI_BCAST(neut_delta_v, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SHARMEM, ierr_mpi)
    CALL MPI_BCAST(neut_delta_w, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SHARMEM, ierr_mpi)
    CALL MPI_BCAST(neut_dir_u, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SHARMEM, ierr_mpi)
    CALL MPI_BCAST(neut_dir_v, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SHARMEM, ierr_mpi)
    CALL MPI_BCAST(neut_dir_w, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SHARMEM, ierr_mpi)
    if (myid_sharmem /= master) ALLOCATE(neut_density(n_u,n_v,n_w),neut_grid_u(n_u),neut_grid_v(n_v),neut_grid_w(n_w))
    CALL MPI_BCAST(neut_grid_u, n_u, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SHARMEM, ierr_mpi)
    CALL MPI_BCAST(neut_grid_v, n_v, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SHARMEM, ierr_mpi)
    CALL MPI_BCAST(neut_grid_w, n_w, MPI_DOUBLE_PRECISION, 0, MPI_COMM_SHARMEM, ierr_mpi)

    CALL mpialloc(neutdens_fspl, 8, n_u, n_v, n_w, myid_sharmem, 0, MPI_COMM_SHARMEM, win_neutdens_spl)

    IF (myid_sharmem == master) THEN
        neutdens_fspl = spline_local%fspl
        CALL EZspline_free(spline_local,ier)
    END IF

    CALL MPI_BARRIER(MPI_COMM_SHARMEM, ier)
#else
    IF (myworkid == master) WRITE(6,*) 'ERROR: HDF5 support required for neutral density reading'   
#endif 

    END SUBROUTINE beams3d_read_neutdens

    SUBROUTINE beams3d_get_neutdens(xq,yq,zq,dens)
        !-----------------------------------------------------------------------
        !     Subroutine: beams3d_get_neutdens
        !           Interpolates the neutdens spline at requested coordinate
        !-----------------------------------------------------------------------
        !     Input Variables
        !          xq,yq,zq     Queried coordinates in x,y,z space
        !     Output variables
        !          dens         Density at position in m^-3
        !-----------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(in) :: xq,yq,zq
        DOUBLE PRECISION, INTENT(out) :: dens
        DOUBLE PRECISION :: pos(3), pos_u, pos_v, pos_w
        INTEGER :: ier

        INTEGER, PARAMETER :: ict(8) = (/1,0,0,0,0,0,0,0/)
        INTEGER :: i,j,k
        DOUBLE PRECISION :: fval(1), xparam, yparam, zparam

        
        ! First remove offset
        pos(1) = xq - neut_x0(1)
        pos(2) = yq - neut_x0(2)
        pos(3) = zq - neut_x0(3)

        ! Get coordinate in u,v,w space
        pos_u = DOT_PRODUCT(pos,neut_dir_u)/neut_delta_u
        pos_v = DOT_PRODUCT(pos,neut_dir_v)/neut_delta_v
        pos_w = DOT_PRODUCT(pos,neut_dir_w)/neut_delta_w
	
        ! Get indices
        i = MAX(INT(pos_u*hri_u),1)
        j = MAX(INT(pos_v*hri_v),1)
        k = MAX(INT(pos_w*hri_w),1)

        ! Get position in [u,v,w]
        xparam = (pos_u - neut_grid_u(i)) * hri_u
        yparam = (pos_v - neut_grid_v(j)) * hri_v
        zparam = (pos_w - neut_grid_w(k)) * hri_w

        CALL R8HERM3FCN(ict, 1, 1, fval, i, j, k, xparam, yparam, zparam, &
                        hr_u, hri_u, hr_v, hri_v, hr_w, hri_w, &
                        neutdens_fspl(1,1,1,1), n_u, n_v, n_w)
        dens = fval(1)
        
    END SUBROUTINE beams3d_get_neutdens


END MODULE beams3d_neutdens
