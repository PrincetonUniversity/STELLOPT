module sources

#if defined MPI && defined ISO_C_BINDING
  use mpi
#endif

  implicit none

  public :: init_sources, finish_sources
  public :: init_quasineutrality_source
  public :: init_source_timeaverage
  public :: update_quasineutrality_source
  public :: include_qn_source
  public :: include_krook_operator, update_tcorr_krook
  public :: remove_zero_projection, project_out_zero
  public :: add_krook_operator
  public :: tcorr_source, exclude_boundary_regions, exp_fac
  public :: tcorr_source_qn, exclude_boundary_regions_qn, exp_fac_qn
  public :: int_krook, int_proj
  public :: qn_source_initialized
#if defined MPI && defined ISO_C_BINDING
  public :: qn_window
#endif

  private

  logical :: include_krook_operator, remove_zero_projection
  logical :: krook_odd, exclude_boundary_regions
  logical :: exclude_boundary_regions_qn
  logical :: from_zero
  real :: nu_krook, tcorr_source, int_krook, int_proj
  real :: tcorr_source_qn
  integer:: ikxmax_source
  real :: exp_fac, exp_fac_qn

  logical :: qn_source_initialized, include_qn_source

#if defined MPI && defined ISO_C_BINDING
  integer :: qn_window = MPI_WIN_NULL
#endif

contains

  subroutine init_sources

    use mp, only: job
    use run_parameters, only: fphi
    use run_parameters, only: ky_solve_radial, ky_solve_real
    use kt_grids, only: nakx, zonal_mode
    use zgrid, only: nzgrid, ntubes
    use stella_layouts, only: vmu_lo
    use dist_fn_arrays, only: g_krook, g_proj
    use fields_arrays, only : phi_proj, phi_proj_stage
    use physics_flags, only: radial_variation
    use species, only: spec, has_electron_species
    use dist_fn, only: adiabatic_option_switch
    use dist_fn, only: adiabatic_option_fieldlineavg
    use file_utils, only: runtype_option_switch, runtype_multibox

    implicit none

    logical :: has_elec, adia_elec
    real :: fac

    call read_parameters

    if (include_krook_operator.and..not.allocated(g_krook)) then
      allocate (g_krook(nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      g_krook = 0.
    endif

    if (remove_zero_projection.and..not.allocated(g_proj)) then
      allocate (g_proj(nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      g_proj = 0.
    endif

    if (.not.allocated(phi_proj)) then
      allocate (phi_proj(nakx,-nzgrid:nzgrid,ntubes)); phi_proj = 0.
    endif
    if (.not.allocated(phi_proj_stage)) then
      allocate (phi_proj_stage(nakx,-nzgrid:nzgrid,ntubes)); phi_proj_stage = 0.
    endif

    fac = 1.
    if (from_zero) fac = 0.

    if (int_krook.lt.0.) int_krook = fac*tcorr_source
    if (int_proj.lt.0.)  int_proj  = fac*tcorr_source

    include_qn_source = .false.
    if (fphi > epsilon(0.0).and.radial_variation.and.ky_solve_radial.gt.0) then
      has_elec  = has_electron_species(spec)
      adia_elec = .not.has_elec.and.zonal_mode(1) &
                  .and.adiabatic_option_switch == adiabatic_option_fieldlineavg
      if(adia_elec) then
        if(runtype_option_switch.ne.runtype_multibox.or.(job.eq.1.and..not.ky_solve_real)) then
          include_qn_source = .true.
        endif
      endif
    endif

  end subroutine init_sources

  subroutine read_parameters

    use file_utils, only: input_unit_exist
    use physics_flags, only: full_flux_surface, radial_variation
    use mp, only: proc0, broadcast
    use kt_grids, only: ikx_max, periodic_variation

    implicit none

    namelist /sources/ &
         include_krook_operator, nu_krook, tcorr_source, remove_zero_projection, &
         ikxmax_source, krook_odd, exclude_boundary_regions, &
         tcorr_source_qn, exclude_boundary_regions_qn, from_zero

    integer :: in_file
    logical :: dexist

    if (proc0) then
       include_krook_operator = .false.
       exclude_boundary_regions = radial_variation.and..not.periodic_variation
       exclude_boundary_regions_qn = exclude_boundary_regions
       remove_zero_projection = .false.
       nu_krook = 0.05
       tcorr_source = 0.02
       tcorr_source_qn = -1.0
       ikxmax_source = 1 ! kx=0
       if(periodic_variation) ikxmax_source = 2 ! kx=0 and kx=1
       krook_odd = .true. ! damp only the odd mode that can affect profiles
       from_zero = .true.

       in_file = input_unit_exist("sources", dexist)
       if (dexist) read (unit=in_file, nml=sources)

       if (tcorr_source_qn.lt.0) tcorr_source_qn = tcorr_source
    end if

    ikxmax_source = min(ikxmax_source,ikx_max)

    int_proj = -1.
    int_krook= -1.

    call broadcast (include_krook_operator)
    call broadcast (exclude_boundary_regions)
    call broadcast (exclude_boundary_regions_qn)
    call broadcast (nu_krook)
    call broadcast (tcorr_source)
    call broadcast (tcorr_source_qn)
    call broadcast (ikxmax_source)
    call broadcast (remove_zero_projection)
    call broadcast (krook_odd)
    call broadcast (from_zero)

  end subroutine read_parameters

  subroutine init_source_timeaverage

    use stella_time, only: code_dt

    implicit none

    exp_fac    = exp(-code_dt/tcorr_source)
    exp_fac_qn = exp(-code_dt/tcorr_source_qn)

  end subroutine init_source_timeaverage

  subroutine finish_sources

    use dist_fn_arrays, only: g_krook, g_proj
    use fields_arrays, only : phi_proj, phi_proj_stage
#if !defined(MPI) || !defined(ISO_C_BINDING)
    use fields_arrays, only : phizf_solve, phi_ext
#endif

    implicit none

    integer :: ierr

    if (allocated(g_krook))  deallocate (g_krook)
    if (allocated(g_proj))   deallocate (g_proj)
    if (allocated(phi_proj)) deallocate (phi_proj)
    if (allocated(phi_proj_stage)) deallocate (phi_proj_stage)

#if defined MPI && defined ISO_C_BINDING
    if (qn_window.ne.MPI_WIN_NULL) call mpi_win_free(qn_window, ierr)
#else    
    if (associated(phizf_solve%zloc)) deallocate (phizf_solve%zloc)
    if (associated(phizf_solve%idx)) deallocate (phizf_solve%idx)
    if (associated(phi_ext)) deallocate (phi_ext)
#endif

  end subroutine finish_sources

  subroutine add_krook_operator (g, gke_rhs)

    use zgrid, only: nzgrid, ntubes
    use constants, only: pi, zi
    use kt_grids, only: akx, nakx, zonal_mode
    use stella_layouts, only: vmu_lo
    use stella_time, only: code_dt
    use dist_fn_arrays, only: g_krook
    use multibox, only: copy_size
    use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded

    implicit none

    complex :: tmp
    integer :: ikx, jkx, iz, it, ia, ivmu, npts

    !complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), optional, intent (in) :: f0
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:),  intent (in) :: g
    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:), intent (in out) :: gke_rhs

    complex, dimension (:,:), allocatable :: g0k, g0x, g1x
    real, dimension (:), allocatable :: basis_func

    ia = 1
    if(.not.zonal_mode(1)) return

    !TODO: add number and momentum conservation
    if (exclude_boundary_regions) then
      npts = nakx - 2*copy_size
      allocate (g0k(1,nakx))
      allocate (g0x(1,nakx))
      allocate (g1x(1,nakx))
      allocate (basis_func(npts))
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
        do it = 1, ntubes
          do iz = -nzgrid, nzgrid
            g0k(1,:) = g(1,:,iz,it,ivmu)
            g1x = 0.
            call transform_kx2x_unpadded(g0k,g0x)
            do ikx = 1, ikxmax_source
              if (ikx.eq.1) then
                basis_func = 1.0
                tmp = sum(g0x(1,(copy_size+1):(nakx-copy_size)))/real(npts)
              else
                do jkx = 1, npts
                  basis_func(jkx) = sin(2.0*pi*(ikx-1)*jkx/real(npts+1))
                enddo
                tmp = 2.0*sum(basis_func*g0x(1,(copy_size+1):(nakx-copy_size)))/real(npts+1)
              endif
              if(tcorr_source.gt.epsilon(0.0)) then
                tmp = (code_dt*tmp + exp_fac*int_krook*g_krook(ikx,iz,it,ivmu)) &
                    / (code_dt     + exp_fac*int_krook)
              endif
              do jkx = 1, npts
                g1x(1,copy_size+jkx) = g1x(1,copy_size+jkx) + tmp*basis_func(jkx) 
              enddo
            enddo
            call transform_x2kx_unpadded (g1x,g0k)
            gke_rhs(1,:,iz,it,ivmu) = gke_rhs(1,:,iz,it,ivmu) - code_dt*nu_krook*g0k(1,:)
          enddo
        enddo
      enddo
      deallocate (g0k, g0x, g1x, basis_func)
    else
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
        do it = 1, ntubes
          do iz = -nzgrid, nzgrid
            do ikx = 1, nakx
              if(abs(akx(ikx)).gt.akx(ikxmax_source)) cycle
              tmp = g(1,ikx,iz,it,ivmu)
              if(krook_odd.and.abs(akx(ikx)).gt.epsilon(0.0)) tmp = zi*aimag(tmp)
              if(tcorr_source.le.epsilon(0.0)) then
                gke_rhs(1,ikx,iz,it,ivmu) = gke_rhs(1,ikx,iz,it,ivmu) - code_dt*nu_krook*tmp
              else
                gke_rhs(1,ikx,iz,it,ivmu) = gke_rhs(1,ikx,iz,it,ivmu) - code_dt*nu_krook &
                                          * (code_dt*tmp + exp_fac*int_krook*g_krook(ikx,iz,it,ivmu)) &
                                          / (code_dt     + exp_fac*int_krook)
              endif
            enddo
          enddo
        enddo
      enddo
    endif

  end subroutine add_krook_operator

  subroutine update_tcorr_krook (g)

    use constants, only: pi, zi
    use dist_fn_arrays, only: g_krook
    use zgrid, only: nzgrid, ntubes
    use kt_grids, only: akx, nakx, zonal_mode
    use stella_layouts, only: vmu_lo
    use stella_time, only: code_dt
    use multibox, only: copy_size
    use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded

    implicit none

    complex, dimension (:,:,-nzgrid:,:,vmu_lo%llim_proc:),  intent (in) :: g
    complex, dimension (:,:), allocatable :: g0k, g0x

    integer :: ivmu, iz, it, ikx, jkx, ia, npts
    real :: int_krook_old
    complex :: tmp

    if(.not.zonal_mode(1)) return

    ia = 1

    int_krook_old = int_krook
    int_krook =  code_dt + exp_fac*int_krook_old

    if (exclude_boundary_regions) then
      npts = nakx - 2*copy_size
      allocate (g0k(1,nakx))
      allocate (g0x(1,nakx))
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
        do it = 1, ntubes
          do iz = -nzgrid, nzgrid
            g0k(1,:) = g(1,:,iz,it,ivmu)
            call transform_kx2x_unpadded(g0k,g0x)
            do ikx = 1, ikxmax_source
              if (ikx.eq.1) then
                tmp = sum(g0x(1,(copy_size+1):(nakx-copy_size)))/real(npts)
              else
                tmp = 0.
                do jkx = 1, npts
                  tmp = tmp + sin(2.0*pi*(ikx-1)*jkx/real(npts+1))*g0x(1,copy_size+jkx)
                enddo
                tmp = 2.0*tmp/real(npts+1)
              endif
            enddo
            g_krook(ikx,iz,it,ivmu) = (code_dt*tmp + exp_fac*int_krook_old*g_krook(ikx,iz,it,ivmu)) &
                                      /int_krook
          enddo
        enddo
      enddo
      deallocate (g0k, g0x)
    else
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
        do it = 1, ntubes
          do iz = -nzgrid, nzgrid
            do ikx = 1, nakx
              tmp = g(1,ikx,iz,it,ivmu)
              if(krook_odd.and.abs(akx(ikx)).gt.epsilon(0.0)) tmp = zi*aimag(tmp)
              g_krook(ikx,iz,it,ivmu) = (code_dt*tmp + exp_fac*int_krook_old*g_krook(ikx,iz,it,ivmu))/int_krook
            enddo
          enddo
        enddo
      enddo
    endif

  end subroutine update_tcorr_krook

  subroutine project_out_zero (g)

    use zgrid, only: nzgrid, ntubes
    use constants, only: pi, zi
    use kt_grids, only: zonal_mode, akx, nakx
    use stella_layouts, only: vmu_lo
    use stella_time, only: code_dt
    use dist_fn_arrays, only: g_proj
    use multibox, only: copy_size
    use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded

    implicit none

    complex :: tmp
    integer :: ikx, jkx, iz, it, ia, ivmu, npts

    complex, dimension (:,:), allocatable :: g0k, g0x, g1x
    real, dimension (:), allocatable :: basis_func
    complex, dimension (:,-nzgrid:,:,vmu_lo%llim_proc:),  intent (inout) :: g

    ia = 1
    if(.not.zonal_mode(1)) return

    if (exclude_boundary_regions) then
      npts = nakx - 2*copy_size
      allocate (g0k(1,nakx))
      allocate (g0x(1,nakx))
      allocate (g1x(1,nakx))
      allocate (basis_func(npts))
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
        do it = 1, ntubes
          do iz = -nzgrid, nzgrid
            g0k(1,:) = g(:,iz,it,ivmu)
            g1x = 0.
            call transform_kx2x_unpadded (g0k,g0x)
            do ikx = 1, ikxmax_source
              !physical region should have an odd number of collocation points
              if (ikx.eq.1) then
                basis_func = 1.0
                tmp = sum(g0x(1,(copy_size+1):(nakx-copy_size)))/real(npts)
              else
                ! here we use a Fourier basis due to periodicity, 
                ! though we could use Legendre polynomials
                ! NB: Only a constant or linear function (or nearly linear, i.e. first
                ! sine harmonic) make physical sense as sources, so ikxmax_source <= 2
                do jkx = 1, npts
                  basis_func(jkx) = sin(2.0*pi*(ikx-1)*jkx/real(npts+1))
                enddo
                tmp = 2.0*sum(basis_func*g0x(1,(copy_size+1):(nakx-copy_size)))/real(npts+1)
              endif
              if(tcorr_source.gt.epsilon(0.)) then
                tmp = (code_dt*tmp + exp_fac*int_proj*g_proj(ikx,iz,it,ivmu)) &
                    / (code_dt     + exp_fac*int_proj)
                g_proj(ikx,iz,it,ivmu) = tmp
              endif
              do jkx = 1, npts
                g1x(1,copy_size+jkx) = g1x(1,copy_size+jkx) + tmp*basis_func(jkx) 
              enddo
            enddo
            call transform_x2kx_unpadded (g1x,g0k)
            g(:,iz,it,ivmu) = g0k(1,:)
          enddo
        enddo
      enddo
      deallocate (g0k, g0x, g1x, basis_func)
    else
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
        do it = 1, ntubes
          do iz = -nzgrid, nzgrid
            do ikx = 1, nakx
              if(abs(akx(ikx)).gt.akx(ikxmax_source)) then
                g(ikx,iz,it,ivmu) = 0.0
              else
                tmp = g(ikx,iz,it,ivmu)
                if(krook_odd.and.abs(akx(ikx)).gt.epsilon(0.0)) tmp = zi*aimag(tmp)
                if(tcorr_source.le.epsilon(0.)) then
                  g(ikx,iz,it,ivmu) = tmp
                else
                  g(ikx,iz,it,ivmu) = (code_dt*tmp + exp_fac*int_proj*g_proj(ikx,iz,it,ivmu)) &
                                    / (code_dt     + exp_fac*int_proj)
                endif
              endif
              if(krook_odd.and.abs(akx(ikx)).gt.epsilon(0.0)) then
                g_proj(ikx,iz,it,ivmu) = zi*aimag(g(ikx,iz,it,ivmu))
              else
                g_proj(ikx,iz,it,ivmu) = g(ikx,iz,it,ivmu)
              endif
            enddo
          enddo
        enddo
      enddo
    endif

    int_proj = code_dt + exp_fac*int_proj

  end subroutine project_out_zero

  subroutine init_quasineutrality_source

    use mp, only: sum_allreduce
#if defined MPI && defined ISO_C_BINDING
    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer
    use mp, only: sgproc0, curr_focus, mp_comm, sharedsubprocs, comm_sgroup
    use mp, only: scope, real_size, nbytes_real
    use mp_lu_decomposition, only: lu_decomposition_local, lu_inverse_local
    use mpi
#endif
    use physics_flags, only: radial_variation
    use stella_geometry, only: dl_over_b, d_dl_over_b_drho
    use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded
    use zgrid, only: nzgrid, nztot
    use kt_grids, only: naky, nakx, rho_d_clamped
    use linear_solve, only: lu_decomposition
    use multibox, only: copy_size
    use fields_arrays, only: phizf_solve, c_mat, theta, phi_ext
    

    implicit none

    integer :: iz, ikx, ia, jkx, jz
    integer :: inmat, jnmat, nmat_zf
    real :: dum
#if defined MPI && ISO_C_BINDING
    integer :: prior_focus, ierr
    integer :: disp_unit = 1
    integer*8 :: cur_pos
    integer (kind=MPI_ADDRESS_KIND) :: win_size
    type(c_ptr) :: cptr
    complex, dimension (:,:), allocatable :: temp_mat
#endif
    complex, dimension (:,:), allocatable :: g0k, g0x, g1k

    ia = 1

    if (qn_source_initialized) return
    qn_source_initialized = .true.

    if (include_qn_source) then
      nmat_zf = nakx*(nztot-1)
#if defined MPI && ISO_C_BINDING
      if(qn_window.eq.MPI_WIN_NULL) then
        prior_focus = curr_focus
        call scope (sharedsubprocs)
        win_size = 0
        if(sgproc0) then
          win_size = int(nmat_zf,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND & 
                   + int(nmat_zf*(nmat_zf+1),MPI_ADDRESS_KIND)*2*real_size !complex size
        endif

        call mpi_win_allocate_shared(win_size,disp_unit,MPI_INFO_NULL, &
                                     mp_comm,cptr,qn_window,ierr)

        if (.not.sgproc0) then
          !make sure all the procs have the right memory address
          call mpi_win_shared_query(qn_window,0,win_size,disp_unit,cptr,ierr)
        end if
        call mpi_win_fence(0,qn_window,ierr)
        cur_pos = transfer(cptr,cur_pos)

        !allocate the memory
        if (.not.associated(phizf_solve%zloc))  then
          cptr = transfer(cur_pos,cptr)
          call c_f_pointer (cptr,phizf_solve%zloc,(/nmat_zf,nmat_zf/))
        endif
        cur_pos = cur_pos + nmat_zf**2*2*nbytes_real

        if (.not.associated(phi_ext))  then
          cptr = transfer(cur_pos,cptr)
          call c_f_pointer (cptr,phi_ext,(/nmat_zf/))
        endif
        cur_pos = cur_pos + nmat_zf*2*nbytes_real

        if (.not.associated(phizf_solve%idx))  then
          cptr = transfer(cur_pos,cptr)
          call c_f_pointer (cptr,phizf_solve%idx,(/nmat_zf/))
        endif

        call mpi_win_fence(0,qn_window,ierr)

        call scope (prior_focus)
      endif
#else
      if (.not.associated(phizf_solve%zloc)) allocate (phizf_solve%zloc(nmat_zf,nmat_zf))
      if (.not.associated(phizf_solve%idx))  allocate (phizf_solve%idx(nmat_zf))
      if (.not.associated(phi_ext)) allocate (phi_ext)
#endif

#if defined MPI && ISO_C_BINDING
      if(sgproc0) then
#endif
        allocate (g0k(1,nakx))
        allocate (g0x(1,nakx))
        allocate (g1k(1,nakx))

        phizf_solve%zloc = 0.

        !get the big matrix
        do jz = -nzgrid, nzgrid-1
          do jkx = 1, nakx
            jnmat = jkx + nakx*(jz+nzgrid)
            ! C.phi
            do ikx = 1, nakx
              inmat = ikx + nakx*(jz+nzgrid)
              phizf_solve%zloc(inmat,jnmat) = phizf_solve%zloc(inmat,jnmat) + c_mat(ikx,jkx)
            enddo

            ! -C.<phi>_\psi
            g0k = 0.0; g0k(1,jkx) = 1.0
            call transform_kx2x_unpadded (g0k,g0x)
            g0x(1,:) = (dl_over_b(ia,jz) + d_dl_over_b_drho(ia,jz)*rho_d_clamped)*g0x(1,:)
            call transform_x2kx_unpadded(g0x,g0k)

            !set the gauge potential
            if(jkx.eq.1) g0k(1,1) = 0. 

            do ikx = 1, nakx
              g1k(1,ikx) = sum(c_mat(ikx,:)*g0k(1,:))
            enddo

            do iz = -nzgrid, nzgrid-1
              do ikx = 1, nakx
                inmat = ikx + nakx*(iz+nzgrid)
                phizf_solve%zloc(inmat,jnmat) = phizf_solve%zloc(inmat,jnmat) - g1k(1,ikx)
              enddo
            enddo

            ! get theta.phi
            g1k(1,:) = theta(:,jkx,jz)

            ! +theta.phi
            do ikx = 1, nakx
              inmat = ikx + nakx*(jz+nzgrid)
              phizf_solve%zloc(inmat,jnmat) = phizf_solve%zloc(inmat,jnmat) + g1k(1,ikx)
            enddo

            ! -<<theta.phi>_psi>_T
            call transform_kx2x_unpadded (g1k, g0x)
            g0x(1,:) = (dl_over_b(ia,jz) + d_dl_over_b_drho(ia,jz)*rho_d_clamped)*g0x(1,:)

            if (exclude_boundary_regions_qn) then
              g0x(1,:) = sum(g0x(1,(copy_size+1):(nakx-copy_size))) &
                       / (nakx - 2*copy_size)
              g0x(1,1:copy_size) = 0.0
              g0x(1,(nakx-copy_size+1):) = 0.0
            else
              g0x(1,:) = sum(g0x(1,:))/nakx
            endif

            call transform_x2kx_unpadded(g0x,g0k)

            if (tcorr_source_qn.gt.epsilon(0.)) then
              g0k = (1. - exp_fac_qn)*g0k
            endif

            do iz = -nzgrid, nzgrid-1
              do ikx = 1, nakx
                inmat = ikx + nakx*(iz+nzgrid)
                phizf_solve%zloc(inmat,jnmat) = phizf_solve%zloc(inmat,jnmat) &
                                              - g0k(1,ikx)
              enddo
            enddo
          enddo
        enddo
        deallocate (g0k,g1k,g0x)
#if defined MPI && ISO_C_BINDING
      endif
      call mpi_win_fence(0,qn_window,ierr)
      call lu_decomposition_local(comm_sgroup, 0, qn_window, &
                                  phizf_solve%zloc, phizf_solve%idx, dum)

      allocate (temp_mat(nmat_zf,nmat_zf))
      temp_mat = phizf_solve%zloc

      call mpi_win_fence(0,qn_window,ierr)

      ! inverse is calculated since it is more straightforward to parallelize
      ! inverse calculation/matrix multiplication than the lu back substitution
      call lu_inverse_local(comm_sgroup, 0, qn_window, &
                            temp_mat, phizf_solve%idx, phizf_solve%zloc)
      deallocate (temp_mat)
#else
      call lu_decomposition(phizf_solve%zloc, phizf_solve%idx, dum)
#endif
    endif
  end subroutine init_quasineutrality_source

  subroutine update_quasineutrality_source

    use fields_arrays, only: phi_proj, phi_proj_stage

    implicit none

    if (tcorr_source_qn.lt.epsilon(0.)) then
      phi_proj = phi_proj_stage
    else
      phi_proj = exp_fac_qn*phi_proj + (1.-exp_fac_qn)*phi_proj_stage
    endif

  end subroutine update_quasineutrality_source

end module sources
