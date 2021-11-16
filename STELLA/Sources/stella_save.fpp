# include "define.inc"

module stella_save

  use mp, only: mp_comm, mp_info

# ifdef NETCDF
!  use netcdf, only: NF90_FLOAT, NF90_DOUBLE
# ifdef NETCDF_PARALLEL
! If using netcdf version 4.1.2 or older delete NF90_MPIIO
  use netcdf, only: NF90_HDF5,NF90_MPIIO
  use netcdf, only: nf90_var_par_access, NF90_COLLECTIVE
  use netcdf, only: nf90_put_att, NF90_GLOBAL, nf90_get_att
# endif
  use netcdf, only: NF90_NOWRITE, NF90_CLOBBER, NF90_NOERR
  use netcdf, only: nf90_create, nf90_open, nf90_sync, nf90_close
  use netcdf, only: nf90_def_dim, nf90_def_var, nf90_enddef
  use netcdf, only: nf90_put_var, nf90_get_var, nf90_strerror
  use netcdf, only: nf90_inq_dimid, nf90_inquire_dimension
  use netcdf, only: nf90_inq_varid, nf90_inquire_variable
  use netcdf, only: nf90_int
  
  use netcdf_utils, only: get_netcdf_code_precision
  use netcdf_utils, only: check_netcdf_file_precision
  use netcdf_utils, only: netcdf_error
  use netcdf_utils, only: netcdf_real, kind_nf
# endif

  implicit none

  public :: stella_restore, stella_save_for_restart
  public :: read_many, save_many
  public :: init_save, init_dt, init_tstart, finish_save

!# ifdef NETCDF
!  public :: netcdf_real, kind_nf, get_netcdf_code_precision, netcdf_error
!# endif

  interface stella_restore
     module procedure stella_restore_many
  end interface

  logical :: read_many = .true., save_many = .true. ! Read and write single or multiple restart files
  
  private
  character (300), save :: restart_file

# ifdef NETCDF
  real, allocatable, dimension (:,:,:) :: tmpr, tmpi
  real, allocatable, dimension (:,:,:,:) :: ktmpr, ktmpi
  real, allocatable, dimension (:,:,:,:)   :: ptmpr, ptmpi
  real, allocatable, dimension (:,:,:)   :: pptmpr, pptmpi
  integer (kind_nf) :: ncid, zedid, vpaid, gloid, gvmuloid, kyid, kxid, muid, tubeid
  integer (kind_nf) :: krookr_id, krooki_id, projr_id, proji_id
  integer (kind_nf) :: phiprojr_id, phiproji_id
!  integer (kind_nf) :: bparr_id, bpari_id
  integer (kind_nf) :: t0id, gr_id, gi_id, delt0id, istep0id
  integer (kind_nf) :: intkrook_id, intproj_id;
  integer (kind_nf) :: shift_id

  logical :: initialized = .false.
# endif

contains

!!----------------------------------------------------------------------!!
!!----------------------------------------------------------------------!!
!!--Save----------------------------------------------------------------!!
!!----------------------------------------------------------------------!!
!!----------------------------------------------------------------------!!

  subroutine stella_save_for_restart &
       (g, istep0, t0, delt0, istatus, exit_in, fileopt)

# ifdef NETCDF
    use fields_arrays, only: shift_state, phi_proj
    use dist_fn_arrays, only: g_krook, g_proj
    use kt_grids, only: naky, nakx
# else
    use mp, only: proc0
# endif    
    use mp, only: iproc, barrier
# ifdef NETCDF_PARALLEL
    use zgrid, only: nztot
# endif
    use zgrid, only: nzgrid, ntubes
    ! Must include kxkyz_layout_type here to avoid obscure bomb while compiling
    ! stella_diagnostics.f90 (which uses this module) with the Compaq F90 compiler:
    use stella_layouts, only: kxkyz_lo, vmu_lo
    use common_types, only: kxkyz_layout_type
    use file_utils, only: error_unit
    use vpamu_grids, only: nvpa, nmu
    use sources, only: include_krook_operator, int_krook
    use sources, only: remove_zero_projection, int_proj
    use sources, only: include_qn_source
    use physics_flags, only: prp_shear_enabled

    implicit none

    complex, dimension (:,:,kxkyz_lo%llim_proc:), intent (in) :: g
    real, intent (in) :: t0, delt0
    integer, intent (in) :: istep0
    integer, intent (out) :: istatus
    logical, intent (in), optional :: exit_in
    character (20), INTENT (in), optional :: fileopt
# ifdef NETCDF
    character (306) :: file_proc
    character (10) :: suffix
    integer :: i, n_elements, nvmulo_elements, ierr
    integer :: total_elements, total_vmulo_elements
    logical :: has_vmulo
# ifdef NETCDF_PARALLEL
    integer, dimension(3) :: start_pos, counts
# endif
    logical :: exit

!*********-----------------------_**********************

    istatus = 0
    if (present(exit_in)) then
       exit = exit_in
    else
       exit = .false.
    end if

!    if (proc0) then
!      write (*,*) "Starting save_for_restart in ", restart_file
!      write (*,*) "List restart files"
!      call system("echo 'start' >> filelist.txt; ls nc/* >> filelist.txt;  ")
!    end if

    n_elements = kxkyz_lo%ulim_proc-kxkyz_lo%llim_proc+1
    total_elements = kxkyz_lo%ulim_world+1

    nvmulo_elements = vmu_lo%ulim_proc-vmu_lo%llim_proc+1
    total_vmulo_elements = vmu_lo%ulim_world+1

    if (n_elements <= 0) return

    has_vmulo = nvmulo_elements.gt.0.or..not.save_many

    if (.not.initialized) then

       initialized = .true.
       
       file_proc = trim(restart_file)
       
!CMR, 5/4/2011: Add optional piece of filename
       IF (PRESENT(fileopt)) THEN
          file_proc=trim(file_proc)//trim(fileopt)
       END IF
!CMRend 

!</HL>  The NETCDF_PARALLEL directives include code for parallel 
!       netcdf using HDF5 to write the output to a single restart file
!       The read_many flag allows the old style multiple file output
# ifdef NETCDF_PARALLEL
       if(save_many) then
# endif
          WRITE (suffix,'(a1,i0)') '.', iproc
# ifdef NETCDF_PARALLEL
       else
          WRITE (suffix,*) ''
       endif
# endif

       file_proc = trim(trim(file_proc)//adjustl(suffix))

# ifdef NETCDF_PARALLEL       
       if(save_many) then
# endif
          istatus = nf90_create (file_proc, NF90_CLOBBER, ncid)
# ifdef NETCDF_PARALLEL
       else
          call barrier
          
          if(iproc .eq. 0) then
             open(unit=tmpunit, file=file_proc)
             close(unit=tmpunit, status='delete')
          end if

          call barrier
! If using netcdf version 4.1.2 or older replace NF90_MPIIO with NF90_CLOBBER
          istatus = nf90_create (file_proc, IOR(NF90_HDF5,NF90_MPIIO), ncid, comm=mp_comm, info=mp_info)
       end if
# endif

       if (istatus /= NF90_NOERR) then
          ierr = error_unit()
          write(ierr,*) "nf90_create error: ", nf90_strerror(istatus)
          goto 1
       end if

# ifdef NETCDF_PARALLEL
       if(.not.save_many) then
          istatus = nf90_put_att(ncid, NF90_GLOBAL, 'xyzs_layout', xyzs_layout)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_put_attr error: ", nf90_strerror(istatus)
             goto 1
          end if
          istatus = nf90_put_att(ncid, NF90_GLOBAL, 'vms_layout', vms_layout)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_put_attr error: ", nf90_strerror(istatus)
             goto 1
          end if
       endif
# endif
       
       if (n_elements > 0) then
          istatus = nf90_def_dim (ncid, "tube", ntubes, tubeid)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_def_dim zed error: ", nf90_strerror(istatus)
             goto 1
          end if

          istatus = nf90_def_dim (ncid, "zed", 2*nzgrid+1, zedid)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_def_dim zed error: ", nf90_strerror(istatus)
             goto 1
          end if

          istatus = nf90_def_dim (ncid, "vpa", nvpa, vpaid)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_def_dim vpa error: ", nf90_strerror(istatus)
             goto 1
          end if

          istatus = nf90_def_dim (ncid, "mu", nmu, muid)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_def_dim mu error: ", nf90_strerror(istatus)
             goto 1
          end if
   
# ifdef NETCDF_PARALLEL                              
          if(save_many) then
# endif
             istatus = nf90_def_dim (ncid, "glo", n_elements, gloid)
# ifdef NETCDF_PARALLEL                    
          else        
             istatus = nf90_def_dim (ncid, "glo", total_elements, gloid)
          endif
# endif
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_def_dim glo error: ", nf90_strerror(istatus)
             goto 1
          end if

# ifdef NETCDF_PARALLEL                              
          if(save_many) then
# endif
            if(nvmulo_elements.gt.0) then
              istatus = nf90_def_dim (ncid, "gvmulo", nvmulo_elements, gvmuloid)
              if (istatus /= NF90_NOERR) then
                ierr = error_unit()
                write(ierr,*) "nf90_def_dim gvmulo error: ", nf90_strerror(istatus)
                goto 1
              endif 
            endif
# ifdef NETCDF_PARALLEL                    
          else        
            istatus = nf90_def_dim (ncid, "gvmulo", total_vmulo_elements, gvmuloid)
            if (istatus /= NF90_NOERR) then
              ierr = error_unit()
              write(ierr,*) "nf90_def_dim gvmulo error: ", nf90_strerror(istatus)
              goto 1
            endif
          endif
# endif
          
          istatus = nf90_def_dim (ncid, "aky", naky, kyid)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_def_dim aky error: ", nf90_strerror(istatus)
             goto 1
          end if
          
          istatus = nf90_def_dim (ncid, "akx", nakx, kxid)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_def_dim akx error: ", nf90_strerror(istatus)
             goto 1
          end if
       end if
       
       if (netcdf_real == 0) netcdf_real = get_netcdf_code_precision()

       istatus = nf90_def_var (ncid, "t0", netcdf_real, t0id)
       if (istatus /= NF90_NOERR) then
          ierr = error_unit()
          write(ierr,*) "nf90_def_var t0 error: ", nf90_strerror(istatus)
          goto 1
       end if

       istatus = nf90_def_var (ncid, "istep0", nf90_int, istep0id)
       if (istatus /= NF90_NOERR) then
          ierr = error_unit()
          write(ierr,*) "nf90_def_var istep0 error: ", nf90_strerror(istatus)
          goto 1
       end if
       
       istatus = nf90_def_var (ncid, "delt0", netcdf_real, delt0id)
       if (istatus /= NF90_NOERR) then
          ierr = error_unit()
          write(ierr,*) "nf90_def_var delt0 error: ", nf90_strerror(istatus)
          goto 1
       end if
       
       if (n_elements > 0) then
          istatus = nf90_def_var (ncid, "gr", netcdf_real, &
               (/ vpaid, muid, gloid /), gr_id)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_def_var g error: ", nf90_strerror(istatus)
             goto 1
          end if
          
          istatus = nf90_def_var (ncid, "gi", netcdf_real, &
               (/ vpaid, muid, gloid /), gi_id)
          if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write(ierr,*) "nf90_def_var g error: ", nf90_strerror(istatus)
             goto 1
          end if
             
          if (include_krook_operator.and.has_vmulo) then
             istatus = nf90_def_var (ncid, "intkrook", netcdf_real, intkrook_id)
             if (istatus /= NF90_NOERR) then
                ierr = error_unit()
                write(ierr,*) "nf90_def_var intkrook error: ", nf90_strerror(istatus)
                goto 1
             end if

             istatus = nf90_def_var (ncid, "krookr", netcdf_real, &
                  (/ kxid, zedid, tubeid, gvmuloid /), krookr_id)
             if (istatus /= NF90_NOERR) then
                ierr = error_unit()
                write(ierr,*) "nf90_def_var apar error: ", nf90_strerror(istatus)
                goto 1
             end if
             
             istatus = nf90_def_var (ncid, "krooki", netcdf_real, &
                  (/ kxid, zedid, tubeid, gvmuloid /), krooki_id)
             if (istatus /= NF90_NOERR) then
                ierr = error_unit()
                write(ierr,*) "nf90_def_var krooki error: ", nf90_strerror(istatus)
                goto 1
             end if

          end if

          if (include_qn_source.and.iproc.eq.0) then
             istatus = nf90_def_var (ncid, "phiprojr", netcdf_real, &
                  (/ kxid, zedid, tubeid /), phiprojr_id)
             if (istatus /= NF90_NOERR) then
                ierr = error_unit()
                write(ierr,*) "nf90_def_var phiprojr error: ", nf90_strerror(istatus)
                goto 1
             end if
             
             istatus = nf90_def_var (ncid, "phiproji", netcdf_real, &
                  (/ kxid, zedid, tubeid /), phiproji_id)
             if (istatus /= NF90_NOERR) then
                ierr = error_unit()
                write(ierr,*) "nf90_def_var phiproji error: ", nf90_strerror(istatus)
                goto 1
             end if
          endif

          if (remove_zero_projection.and.has_vmulo) then
             istatus = nf90_def_var (ncid, "intproj", netcdf_real, intproj_id)
             if (istatus /= NF90_NOERR) then
                ierr = error_unit()
                write(ierr,*) "nf90_def_var intproj error: ", nf90_strerror(istatus)
                goto 1
             end if

             istatus = nf90_def_var (ncid, "projr", netcdf_real, &
                  (/ kxid, zedid, tubeid, gvmuloid /), projr_id)
             if (istatus /= NF90_NOERR) then
                ierr = error_unit()
                write(ierr,*) "nf90_def_var projr error: ", nf90_strerror(istatus)
                goto 1
             end if
             
             istatus = nf90_def_var (ncid, "proji", netcdf_real, &
                  (/ kxid, zedid, tubeid, gvmuloid /), proji_id)
             if (istatus /= NF90_NOERR) then
                ierr = error_unit()
                write(ierr,*) "nf90_def_var proji error: ", nf90_strerror(istatus)
                goto 1
             end if

          end if

          if (prp_shear_enabled) then
             istatus = nf90_def_var (ncid, "shiftstate", netcdf_real,&
                (/ kyid /), shift_id)
             if (istatus /= NF90_NOERR) then
                ierr = error_unit()
                write(ierr,*) "nf90_def_var shiftstate error: ", nf90_strerror(istatus)
                goto 1
             end if
          endif

!           if (fbpar > epsilon(0.)) then
!              istatus = nf90_def_var (ncid, "bpar_r", netcdf_real, &
!                   (/ zedid, kxid, kyid /), bparr_id)
!              if (istatus /= NF90_NOERR) then
!                 ierr = error_unit()
!                 write(ierr,*) "nf90_def_var bparr error: ", nf90_strerror(istatus)
!                 goto 1
!              end if

!              istatus = nf90_def_var (ncid, "bpar_i", netcdf_real, &
!                   (/ zedid, kxid, kyid /), bpari_id)
!              if (istatus /= NF90_NOERR) then
!                 ierr = error_unit()
!                 write(ierr,*) "nf90_def_var bpari error: ", nf90_strerror(istatus)
!                 goto 1
!              end if
!           end if
          
       end if

! remove allocated conditional because we want to be able to restart
! using exb shear from a case which does not have exb shear (i.e.
! we need kx_shift variable defined in netcdf file even if no exb
! shear present in simulation) -- MAB + CMR
!       if (allocated(kx_shift)) then   ! MR begin
!       istatus = nf90_def_var (ncid, "kx_shift", netcdf_real, &
!            (/ kyid /), kx_shift_id)
!       if (istatus /= NF90_NOERR) then
!          ierr = error_unit()
!          write(ierr,*) "nf90_def_var kx_shift error: ", nf90_strerror(istatus)
!          goto 1
!       endif
!       endif   ! MR end 
        
!    if (proc0) then
!      write (*,*) "Finished definitions"
    !      write (*,*) "List restart files"
    !      call system("echo 'defs' >> filelist.txt; ls nc/* >> filelist.txt;  ")
!    end if
       
       istatus = nf90_enddef (ncid)
       if (istatus /= NF90_NOERR) then
          ierr = error_unit()
          write (ierr,*) "nf90_enddef error: ", nf90_strerror(istatus)
          goto 1
       end if
    end if


    !!!-----------------------!!!
    !!!-----------------------!!!
    !!!-----------------------!!!

# ifdef NETCDF_PARALLEL                    
    if(save_many .or. iproc == 0) then
# endif

       istatus = nf90_put_var (ncid, delt0id, delt0)
       if (istatus /= NF90_NOERR) then
          ierr = error_unit()
          write (ierr,*) "nf90_put_var delt0 error: ", nf90_strerror(istatus)
          goto 1
       end if
 
       istatus = nf90_put_var (ncid, t0id, t0)
       if (istatus /= NF90_NOERR) then
          ierr = error_unit()
          write (ierr,*) "nf90_put_var t0 error: ", nf90_strerror(istatus)
          goto 1
       end if

       istatus = nf90_put_var (ncid, istep0id, istep0)
       if (istatus /= NF90_NOERR) then
          ierr = error_unit()
          write (ierr,*) "nf90_put_var istep0 error: ", nf90_strerror(istatus)
          goto 1
       end if

# ifdef NETCDF_PARALLEL
    endif
# endif

1   continue

    if (istatus /= NF90_NOERR) then
       i = nf90_close (ncid)
       return
    end if

    if (n_elements > 0) then

       if (.not. allocated(tmpr)) &
            allocate (tmpr(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
       
       tmpr = real(g)

# ifdef NETCDF_PARALLEL
       if(save_many) then
# endif
          istatus = nf90_put_var (ncid, gr_id, tmpr)
#ifdef NETCDF_PARALLEL
       else
          istatus = nf90_var_par_access(ncid, gr_id, NF90_COLLECTIVE)
          istatus = nf90_var_par_access(ncid, gi_id, NF90_COLLECTIVE)

          start_pos = (/1,1,kxkyz_lo%llim_proc+1/)
          counts = (/nvpa, nmu, n_elements/)

          istatus = nf90_put_var (ncid, gr_id, tmpr, start=start_pos, count=counts)
       endif
# endif     
       
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, gr_id)
       
       tmpr = aimag(g)
# ifdef NETCDF_PARALLEL
       if(save_many) then
# endif
          istatus = nf90_put_var (ncid, gi_id, tmpr)
#ifdef NETCDF_PARALLEL
       else
          istatus = nf90_put_var (ncid, gi_id, tmpr, start=start_pos, count=counts)
       endif
# endif     

       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, gi_id)
          
       if (include_krook_operator) then
         if (.not. allocated(ktmpr)) &
           allocate (ktmpr(nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         if (.not. allocated(ktmpi)) &
           allocate (ktmpi(nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))

# ifdef NETCDF_PARALLEL                    
         if(save_many .or. iproc == 0) then
# endif

           istatus = nf90_put_var (ncid, intkrook_id, int_krook)
           if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write (ierr,*) "nf90_put_var int_krook error: ", nf90_strerror(istatus)
             goto 1
           end if

# ifdef NETCDF_PARALLEL
         endif
# endif

         ktmpr = real(g_krook)
         ktmpi = aimag(g_krook)
         
# ifdef NETCDF_PARALLEL
         if(save_many) then
# endif
           istatus = nf90_put_var (ncid, krookr_id, ktmpr)
#ifdef NETCDF_PARALLEL
         else
           istatus = nf90_var_par_access(ncid, krookr_id, NF90_COLLECTIVE)
           istatus = nf90_var_par_access(ncid, krooki_id, NF90_COLLECTIVE)

           start_pos = (/1,1,1,vmu_lo%llim_proc+1/)
           counts = (/nakx, nztot, ntubes, nvmulo_elements/)

           istatus = nf90_put_var (ncid, krookr_id, ktmpr, start=start_pos, count=counts)
         endif
# endif     
         if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, krookr_id)



# ifdef NETCDF_PARALLEL
         if(save_many) then
# endif
           istatus = nf90_put_var (ncid, krooki_id, ktmpi)
#ifdef NETCDF_PARALLEL
         else
           istatus = nf90_put_var (ncid, krooki_id, ktmpi, start=start_pos, count=counts)
         endif
# endif     
         if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, krooki_id)
       
       end if

       if (remove_zero_projection) then
         if (.not. allocated(ptmpr)) &
           allocate (ptmpr(nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         if (.not. allocated(ptmpi)) &
           allocate (ptmpi(nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))

# ifdef NETCDF_PARALLEL                    
         if(save_many .or. iproc == 0) then
# endif

           istatus = nf90_put_var (ncid, intproj_id, int_proj)
           if (istatus /= NF90_NOERR) then
             ierr = error_unit()
             write (ierr,*) "nf90_put_var int_proj error: ", nf90_strerror(istatus)
             goto 1
           end if

# ifdef NETCDF_PARALLEL
         endif
# endif

         ptmpr = real(g_proj)
         ptmpi = aimag(g_proj)
         
# ifdef NETCDF_PARALLEL
         if(save_many) then
# endif
           istatus = nf90_put_var (ncid, projr_id, ptmpr)
#ifdef NETCDF_PARALLEL
         else
           istatus = nf90_var_par_access(ncid, projr_id, NF90_COLLECTIVE)
           istatus = nf90_var_par_access(ncid, proji_id, NF90_COLLECTIVE)

           start_pos = (/1,1,1,vmu_lo%llim_proc+1/)
           counts = (/nakx, nztot, ntubes, nvmulo_elements/)

           istatus = nf90_put_var (ncid, projr_id, ptmpr, start=start_pos, count=counts)
         endif
# endif     
         if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, projr_id)

# ifdef NETCDF_PARALLEL
         if(save_many) then
# endif
           istatus = nf90_put_var (ncid, proji_id, ptmpi)
#ifdef NETCDF_PARALLEL
         else
           istatus = nf90_put_var (ncid, proji_id, ptmpi, start=start_pos, count=counts)
         endif
# endif     
         if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, proji_id)
       
       end if

       if (prp_shear_enabled) then
         istatus = nf90_put_var (ncid, shift_id, shift_state)
         if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, shift_id)
       end if

       if (include_qn_source.and.iproc.eq.0) then
         if (.not. allocated(pptmpr)) &
           allocate (pptmpr(nakx,-nzgrid:nzgrid,ntubes))
         if (.not. allocated(pptmpi)) &
           allocate (pptmpi(nakx,-nzgrid:nzgrid,ntubes))

         pptmpr = real(phi_proj)
         pptmpi = aimag(phi_proj)
         
# ifdef NETCDF_PARALLEL
         if(save_many) then
# endif
           istatus = nf90_put_var (ncid, phiprojr_id, pptmpr)
#ifdef NETCDF_PARALLEL
         else
           istatus = nf90_var_par_access(ncid, phiprojr_id, NF90_COLLECTIVE)
           istatus = nf90_var_par_access(ncid, phiproji_id, NF90_COLLECTIVE)

           start_pos = (/1,1,1/)
           counts = (/nakx, nztot, ntubes/)

           istatus = nf90_put_var (ncid, phiprojr_id, ktmpr, start=start_pos, count=counts)
         endif
# endif     
         if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, phiprojr_id)



# ifdef NETCDF_PARALLEL
         if(save_many) then
# endif
           istatus = nf90_put_var (ncid, phiproji_id, pptmpi)
#ifdef NETCDF_PARALLEL
         else
           istatus = nf90_put_var (ncid, phiproji_id, pptmpi, start=start_pos, count=counts)
         endif
# endif     
         if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, phiproji_id)
       
       end if
    end if
       
    if (exit) then
       i = nf90_close (ncid)
       if (i /= NF90_NOERR) &
            call netcdf_error (istatus, message='nf90_close error')
    else
       i = nf90_sync (ncid)
       if (i /= NF90_NOERR) &
            call netcdf_error (istatus, message='nf90_sync error')
    end if

# else

    if (proc0) write (error_unit(),*) &
         'WARNING: stella_save_for_restart is called without netcdf library'

# endif

    if (allocated(tmpr))  deallocate (tmpr)
    if (allocated(tmpi))  deallocate (tmpi)
    if (allocated(ptmpr)) deallocate (ptmpr)
    if (allocated(ptmpi)) deallocate (ptmpi)
    if (allocated(ktmpr)) deallocate (ktmpr)
    if (allocated(ktmpi)) deallocate (ktmpi)
    if (allocated(pptmpr)) deallocate (pptmpr)
    if (allocated(pptmpi)) deallocate (pptmpi)

  end subroutine stella_save_for_restart



!!----------------------------------------------------------------------!!
!!----------------------------------------------------------------------!!
!!---Restart------------------------------------------------------------!!
!!----------------------------------------------------------------------!!
!!----------------------------------------------------------------------!!

  subroutine stella_restore_many (g, scale, istatus)
# ifdef NETCDF
    use fields_arrays, only: shift_state, phi_proj
    use dist_fn_arrays, only: g_krook, g_proj
    use kt_grids, only: naky, nakx
# endif
# ifdef NETCDF_PARALLEL
    use zgrid, only: nztot
# endif
    use mp, only: iproc, broadcast
    use zgrid, only: nzgrid, ntubes
    use vpamu_grids, only: nvpa, nmu
    use stella_layouts, only: kxkyz_lo, vmu_lo
    use file_utils, only: error_unit
    use sources, only: include_krook_operator, int_krook
    use sources, only: remove_zero_projection, int_proj
    use physics_flags, only: prp_shear_enabled
    use sources, only: include_qn_source

    implicit none

    complex, dimension (:,:,kxkyz_lo%llim_proc:), intent (out) :: g
    real, intent (in) :: scale
    integer, intent (out) :: istatus
# ifdef NETCDF
# ifdef NETCDF_PARALLEL
    integer, dimension(3) :: counts, start_pos
# endif
    character (306) :: file_proc
    character (10) :: suffix
    integer :: i, n_elements, nvmulo_elements, ierr
    logical :: has_vmulo
    
    n_elements = kxkyz_lo%ulim_proc-kxkyz_lo%llim_proc+1
    nvmulo_elements = vmu_lo%ulim_proc-vmu_lo%llim_proc+1

    if (n_elements <= 0) return

    has_vmulo = nvmulo_elements.gt.0.or..not.read_many
    
    if (.not.initialized) then
!       initialized = .true.
       file_proc = trim(restart_file)

# ifdef NETCDF_PARALLEL
       if(read_many) then
# endif
          write (suffix,'(a1,i0)') '.', iproc
          file_proc = trim(trim(file_proc)//adjustl(suffix))       
          istatus = nf90_open (file_proc, NF90_NOWRITE, ncid)
# ifdef NETCDF_PARALLEL
       else
! If using netcdf version 4.1.2 deleted NF90_MPIIO and the associated IOR
          istatus = nf90_open (file_proc, IOR(NF90_NOWRITE, NF90_MPIIO), ncid, comm=mp_comm, info=mp_info)
       endif
# endif

       if (istatus /= NF90_NOERR)  then 
         call netcdf_error (istatus, file=file_proc, abort=.true.)
       endif

       ! check precision
       if (netcdf_real == 0) netcdf_real = get_netcdf_code_precision()
       call check_netcdf_file_precision (ncid)

       istatus = nf90_inq_dimid (ncid, "tube", tubeid)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='tube')

       istatus = nf90_inq_dimid (ncid, "zed", zedid)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='zed')
       
       istatus = nf90_inq_dimid (ncid, "aky", kyid)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='aky')
       
       istatus = nf90_inq_dimid (ncid, "akx", kxid)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='akx')
       
       istatus = nf90_inq_dimid (ncid, "glo", gloid)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='glo')

       if(has_vmulo) then
         istatus = nf90_inq_dimid (ncid, "gvmulo", gvmuloid)
         if (istatus /= NF90_NOERR) call netcdf_error (istatus, dim='gvmulo')
       endif
              
       istatus = nf90_inquire_dimension (ncid, tubeid, len=i)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=tubeid)
       if (i /= ntubes) write(*,*) 'Restart error: ntubes=? ',i,' : ',ntubes,' : ',iproc

       istatus = nf90_inquire_dimension (ncid, zedid, len=i)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=zedid)
       if (i /= 2*nzgrid + 1) write(*,*) 'Restart error: nzgrid=? ',i,' : ',nzgrid,' : ',iproc

       istatus = nf90_inquire_dimension (ncid, kyid, len=i)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=kyid)
       if (i /= naky) write(*,*) 'Restart error: naky=? ',i,' : ',naky,' : ',iproc
       
       istatus = nf90_inquire_dimension (ncid, kxid, len=i)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=kxid)
       if (i /= nakx) write(*,*) 'Restart error: nakx=? ',i,' : ',nakx,' : ',iproc
       
       istatus = nf90_inquire_dimension (ncid, gloid, len=i)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=gloid)
#ifdef NETCDF_PARALLEL       
       if(read_many) then
#endif
          if (i /= n_elements) write(*,*) 'Restart error: glo=? ',i,' : ',iproc
#ifdef NETCDF_PARALLEL
       else
          if (i /= kxkyz_lo%ulim_world+1) write(*,*) 'Restart error: glo=? ',i,' : ',iproc
       endif
#endif

       if(has_vmulo) then
         istatus = nf90_inquire_dimension (ncid, gvmuloid, len=i)
         if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, dimid=gvmuloid)
#ifdef NETCDF_PARALLEL       
         if(read_many) then
#endif
           if (i /= nvmulo_elements) write(*,*) 'Restart error: gvmulo=? ',i,' : ',iproc
#ifdef NETCDF_PARALLEL
         else
           if (i /= vmu_lo%ulim_world+1) write(*,*) 'Restart error: gvmulo=? ',i,' : ',iproc
         endif
#endif
       endif
       
       if(include_krook_operator.and.has_vmulo) then
          istatus = nf90_inq_varid (ncid, "intkrook", intkrook_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='intkrook')

          istatus = nf90_inq_varid (ncid, "krookr", krookr_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='krookr')
          
          istatus = nf90_inq_varid (ncid, "krooki", krooki_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='krooki')

       endif

       if(remove_zero_projection.and.has_vmulo) then
          istatus = nf90_inq_varid (ncid, "intproj", intproj_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='intproj')

          istatus = nf90_inq_varid (ncid, "projr", projr_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='projr')
          
          istatus = nf90_inq_varid (ncid, "proji", proji_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='proji')
       endif

       if(include_qn_source.and.iproc.eq.0) then
          istatus = nf90_inq_varid (ncid, "phiprojr", phiprojr_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='phiprojr')
          
          istatus = nf90_inq_varid (ncid, "phiproji", phiproji_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='phiproji')
       endif

       if(prp_shear_enabled) then
          istatus = nf90_inq_varid (ncid, "shiftstate", shift_id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='shiftstate')
       endif

!        if (fbpar > epsilon(0.)) then
!           istatus = nf90_inq_varid (ncid, "bpar_r", bparr_id)
!           if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='bpar_r')

!           istatus = nf90_inq_varid (ncid, "bpar_i", bpari_id)
!           if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='bpar_i')
!        end if

!       if (allocated(kx_shift)) then   ! MR begin
!          istatus = nf90_inq_varid (ncid, "kx_shift", kx_shift_id)
!          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='kx_shift')
!       endif   ! MR end

       istatus = nf90_inq_varid (ncid, "gr", gr_id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='gr')
       
       istatus = nf90_inq_varid (ncid, "gi", gi_id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='gi')
    end if
    
    if (.not. allocated(tmpr)) &
         allocate (tmpr(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
    if (.not. allocated(tmpi)) &
         allocate (tmpi(nvpa,nmu,kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))

    tmpr = 0.; tmpi = 0.
# ifdef NETCDF_PARALLEL
    if(read_many) then
# endif
       istatus = nf90_get_var (ncid, gr_id, tmpr)
#ifdef NETCDF_PARALLEL
    else
       start_pos = (/1,1,kxkyz_lo%llim_proc+1/)
       counts = (/nvpa, nmu, n_elements/)
       istatus = nf90_get_var (ncid, gr_id, tmpr, start=start_pos, count=counts)
    end if
# endif

   if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, gr_id)

# ifdef NETCDF_PARALLEL
    if(read_many) then
# endif
       istatus = nf90_get_var (ncid, gi_id, tmpi)
#ifdef NETCDF_PARALLEL
    else
       istatus = nf90_get_var (ncid, gi_id, tmpi, start=start_pos, count=counts)
    end if
# endif

    if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, gi_id)

    g = cmplx(tmpr, tmpi)
    

    if(include_krook_operator.and.has_vmulo) then
      if (.not. allocated(ktmpr)) &
        allocate (ktmpr(nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      if (.not. allocated(ktmpi)) &
        allocate (ktmpi(nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))

      istatus = nf90_get_var (ncid, intkrook_id, int_krook)
      if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, intkrook_id)

      ktmpr = 0.; ktmpi = 0.
# ifdef NETCDF_PARALLEL
      if(read_many) then
# endif
        istatus = nf90_get_var (ncid, krookr_id, ktmpr)
#ifdef NETCDF_PARALLEL
      else
        start_pos = (/1,1,1,vmu_lo%llim_proc+1/)
        counts = (/nakx, nztot, ntubes, nvmulo_elements/)
        istatus = nf90_get_var (ncid, krookr_id, ktmpr, start=start_pos, count=counts)
      end if
# endif

       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, krookr_id)

# ifdef NETCDF_PARALLEL
      if(read_many) then
# endif
        istatus = nf90_get_var (ncid, krooki_id, ktmpi)
#ifdef NETCDF_PARALLEL
      else
        istatus = nf90_get_var (ncid, krooki_id, ktmpi, start=start_pos, count=counts)
      end if
# endif

      if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, krooki_id)

      g_krook = cmplx(ktmpr, ktmpi)
      
    endif

    if(remove_zero_projection.and.has_vmulo) then
      if (.not. allocated(ptmpr)) &
        allocate (ptmpr(nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      if (.not. allocated(ptmpi)) &
        allocate (ptmpi(nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))

      istatus = nf90_get_var (ncid, intproj_id, int_proj)
      if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, intproj_id)

      ptmpr = 0.; ptmpi = 0.
# ifdef NETCDF_PARALLEL
      if(read_many) then
# endif
        istatus = nf90_get_var (ncid, projr_id, ptmpr)
#ifdef NETCDF_PARALLEL
      else
        start_pos = (/1,1,1,vmu_lo%llim_proc+1/)
        counts = (/nakx, nztot, ntubes, nvmulo_elements/)
        istatus = nf90_get_var (ncid, projr_id, ptmpr, start=start_pos, count=counts)
      end if
# endif

       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, projr_id)

# ifdef NETCDF_PARALLEL
      if(read_many) then
# endif
        istatus = nf90_get_var (ncid, proji_id, ptmpi)
#ifdef NETCDF_PARALLEL
      else
        istatus = nf90_get_var (ncid, proji_id, ptmpi, start=start_pos, count=counts)
      end if
# endif

      if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, proji_id)

      g_proj = cmplx(ptmpr, ptmpi)
      
    endif

    if(include_qn_source.and.iproc.eq.0) then
      if (.not. allocated(pptmpr)) allocate (pptmpr(nakx,-nzgrid:nzgrid,ntubes))
      if (.not. allocated(pptmpi)) allocate (pptmpi(nakx,-nzgrid:nzgrid,ntubes))


      pptmpr = 0.; pptmpi = 0.
# ifdef NETCDF_PARALLEL
      if(read_many) then
# endif
        istatus = nf90_get_var (ncid, phiprojr_id, pptmpr)
#ifdef NETCDF_PARALLEL
      else
        start_pos = (/1,1,1/)
        counts = (/nakx, nztot, ntubes/)
        istatus = nf90_get_var (ncid, phiprojr_id, pptmpr, start=start_pos, count=counts)
      end if
# endif

       if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, phiprojr_id)

# ifdef NETCDF_PARALLEL
      if(read_many) then
# endif
        istatus = nf90_get_var (ncid, phiproji_id, pptmpi)
#ifdef NETCDF_PARALLEL
      else
        istatus = nf90_get_var (ncid, phiproji_id, pptmpi, start=start_pos, count=counts)
      end if
# endif

      if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, phiproji_id)

      phi_proj = cmplx(pptmpr, pptmpi)

    endif

    if(prp_shear_enabled) then
      istatus = nf90_get_var (ncid, shift_id, shift_state)
      if (istatus /= NF90_NOERR) call netcdf_error (istatus, ncid, shift_id)
    endif



    if (scale > 0.) then
       g = g*scale
       if(include_krook_operator) g_krook = g_krook*scale
       if(remove_zero_projection) g_proj = g_proj*scale
    endif

    ! RN 2008/05/23: this was commented out. why? HJL 2013/05/15 Because it stops future writing to the file
!    istatus = nf90_close (ncid)
    if (istatus /= NF90_NOERR) then
       ierr = error_unit()
       write(ierr,*) "nf90_close error: ", nf90_strerror(istatus),' ',iproc
    end if

# else
    
    write (error_unit(),*) &
         'ERROR: stella_restore_many is called without netcdf'

# endif

    if (allocated(tmpr))  deallocate (tmpr)
    if (allocated(tmpi))  deallocate (tmpi)
    if (allocated(ptmpr)) deallocate (ptmpr)
    if (allocated(ptmpi)) deallocate (ptmpi)
    if (allocated(ktmpr)) deallocate (ktmpr)
    if (allocated(ktmpi)) deallocate (ktmpi)
    if (allocated(pptmpr)) deallocate (pptmpr)
    if (allocated(pptmpi)) deallocate (pptmpi)

    if (include_qn_source) call broadcast (phi_proj)

  end subroutine stella_restore_many


  subroutine init_save (file)

    character(300), intent (in) :: file
    
    restart_file = file

  end subroutine init_save

  subroutine init_dt (delt0, istatus)

# ifdef NETCDF
    use mp, only: proc0, broadcast
    use file_utils, only: error_unit
# endif
    implicit none
    real, intent (in out) :: delt0
    integer, intent (out) :: istatus
# ifdef NETCDF
    character (306) :: file_proc

    if (proc0) then

       if (.not. initialized) then

# ifdef NETCDF_PARALLEL
          if(read_many) then
# endif
             file_proc=trim(trim(restart_file)//'.0')
# ifdef NETCDF_PARALLEL
          else 
             file_proc=trim(trim(restart_file))
          end if
# endif

          istatus = nf90_open (file_proc, NF90_NOWRITE, ncid)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus,file=file_proc)

          istatus = nf90_inq_varid (ncid, "delt0", delt0id)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='delt0')
       end if

       istatus = nf90_get_var (ncid, delt0id, delt0)

       if (istatus /= NF90_NOERR) then
          call netcdf_error (istatus, ncid, delt0id, message=' in init_dt')
          delt0 = -1.
       endif           

       if (.not.initialized) istatus = nf90_close (ncid)
    endif

    call broadcast (istatus)
    call broadcast (delt0)

# endif

  end subroutine init_dt

  subroutine init_tstart (tstart, istep0, istatus)

# ifdef NETCDF
    use mp, only: proc0, broadcast
    use file_utils, only: error_unit
# endif
    implicit none
    real, intent (in out) :: tstart
    integer, intent (out) :: istep0
    integer, intent (out) :: istatus
# ifdef NETCDF
    character (306) :: file_proc

    if (proc0) then
# ifdef NETCDF_PARALLEL
       if(read_many) then
# endif
          file_proc=trim(trim(restart_file)//'.0')
# ifdef NETCDF_PARALLEL
       else 
          file_proc=trim(trim(restart_file))
       end if
# endif

       if (.not.initialized) then

          istatus = nf90_open (file_proc, NF90_NOWRITE, ncid)
          if (istatus /= NF90_NOERR) call netcdf_error (istatus, file=file_proc)
       end if

       istatus = nf90_inq_varid (ncid, "t0", t0id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='t0')

       istatus = nf90_get_var (ncid, t0id, tstart)
       if (istatus /= NF90_NOERR) then
          call netcdf_error (istatus, ncid, t0id, message=' in init_tstart')
          tstart = -1.
       end if           

       istatus = nf90_inq_varid (ncid, "istep0", istep0id)
       if (istatus /= NF90_NOERR) call netcdf_error (istatus, var='istep0')

       istatus = nf90_get_var (ncid, istep0id, istep0)
       if (istatus /= NF90_NOERR) then
          call netcdf_error (istatus, ncid, istep0id, message=' in init_tstart')
          istep0 = -1
       end if           

       if (.not.initialized) istatus = nf90_close (ncid)
          
    endif

    call broadcast (istatus)
    call broadcast (istep0)
    call broadcast (tstart)

# endif

  end subroutine init_tstart

  subroutine finish_save
    
    if (allocated(tmpr))  deallocate (tmpr)
    if (allocated(tmpi))  deallocate (tmpi)

  end subroutine finish_save

end module stella_save
