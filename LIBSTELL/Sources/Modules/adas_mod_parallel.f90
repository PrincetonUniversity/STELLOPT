!-----------------------------------------------------------------------
!     Module:        adas_mod_parallel
!     Authors:       S. Lazerson (samuel.lazeson@ipp.mpg.de)
!     Date:          11/21/2019
!     Description:   This module provides a memory efficient, and file
!                    access minimal way to use the ADAS reaction tables. 
!-----------------------------------------------------------------------
module adas_mod_parallel
!-----------------------------------------------------------------------
!     Modules
!-----------------------------------------------------------------------
#if defined(NETCDF)
       USE netcdf
#endif
#if defined(MPI_OPT)
       USE mpi_sharmem
       USE mpi_inc
#endif
!-----------------------------------------------------------------------
!     Variables
!-----------------------------------------------------------------------
       IMPLICIT NONE
       INTEGER :: iwarn_adas
       DOUBLE PRECISION, POINTER, DIMENSION(:) :: ei_1_axis, ei_2_axis
       DOUBLE PRECISION, POINTER, DIMENSION(:) :: ei_1_sigv, ei_2_sigv
       DOUBLE PRECISION, POINTER, DIMENSION(:) :: cx_1_1_axis, cx_2_1_axis,cx_2_2_axis,&
                                                      ii_1_1_axis, ii_1_2_axis, &
                                                      ii_1_3_axis, ii_1_4_axis, &
                                                      ii_1_5_axis, ii_1_6_axis, &
                                                      ii_1_7_axis, ii_1_8_axis, &
                                                      ii_1_9_axis, ii_1_10_axis, &
                                                      ii_2_1_axis, ii_2_2_axis, &
                                                      ii_2_3_axis, ii_2_4_axis, &
                                                      ii_2_5_axis, ii_2_6_axis, &
                                                      ii_2_7_axis, ii_2_8_axis, &
                                                      ex_1_1_e, ex_1_1_n, ex_1_1_t, & 
                                                      ex_1_2_e, ex_1_2_n, ex_1_2_t, &
                                                      ex_1_3_e, ex_1_3_n, ex_1_3_t
       DOUBLE PRECISION, POINTER, DIMENSION(:,:) :: cx_1_1_btsigv,cx_2_1_btsigv,cx_2_2_btsigv,&
                                                        ii_1_1_btsigv, ii_1_2_btsigv, &
                                                        ii_1_3_btsigv, ii_1_4_btsigv, &
                                                        ii_1_5_btsigv, ii_1_6_btsigv, &
                                                        ii_1_7_btsigv, ii_1_8_btsigv, &
                                                        ii_1_9_btsigv, ii_1_10_btsigv, &
                                                        ii_2_1_btsigv, ii_2_2_btsigv, &
                                                        ii_2_3_btsigv, ii_2_4_btsigv, &
                                                        ii_2_5_btsigv, ii_2_6_btsigv, &
                                                        ii_2_7_btsigv, ii_2_8_btsigv
       DOUBLE PRECISION, POINTER, DIMENSION(:,:,:) :: ex_1_1, ex_1_2, ex_1_3
       INTEGER :: win_ei_1a, win_ei_2a, win_cx_11a, win_cx_21a, win_cx_22a,&
                  win_ii_11a, win_ii_12a, win_ii_13a, win_ii_14a, win_ii_15a, &
                  win_ii_16a, win_ii_17a, win_ii_18a, win_ii_19a, win_ii_110a, &
                  win_ii_21a, win_ii_22a, win_ii_23a, win_ii_24a, win_ii_25a, &
                  win_ii_26a, win_ii_27a, win_ii_28a, &
                  win_ei_1s, win_ei_2s, win_cx_11b, win_cx_21b, win_cx_22b,&
                  win_ii_11b, win_ii_12b, win_ii_13b, win_ii_14b, win_ii_15b, &
                  win_ii_16b, win_ii_17b, win_ii_18b, win_ii_19b, win_ii_110b, &
                  win_ii_21b, win_ii_22b, win_ii_23b, win_ii_24b, win_ii_25b, &
                  win_ii_26b, win_ii_27b, win_ii_28b, &
                  win_ex_11_e, win_ex_11_n, win_ex_11_t, &
                  win_ex_12_e, win_ex_12_n, win_ex_12_t, &
                  win_ex_13_e, win_ex_13_n, win_ex_13_t, &
                  win_ex_11, win_ex_12, win_ex_13, win_ex_14, win_ex_15, &
                  win_ex_16, win_ex_17, win_ex_18, win_ex_19, win_ex_110
       INTEGER, PRIVATE, PARAMETER :: local_master=0
!-----------------------------------------------------------------------
!     INTERFACES
!-----------------------------------------------------------------------
       PUBLIC :: adas_btsigv
       PUBLIC :: adas_sigvte_ioniz
       PUBLIC :: adas_bms

       INTERFACE adas_sigvte_ioniz
          MODULE PROCEDURE adas_sigvte_ioniz_r8,adas_sigvte_ioniz_int
       END INTERFACE

       INTERFACE adas_btsigv
          MODULE PROCEDURE adas_btsigv_int_int,adas_btsigv_int_r8,adas_btsigv_r8_r8,adas_btsigv_r8_int
       END INTERFACE

       INTERFACE adas_bms
          MODULE PROCEDURE adas_bms_int_r8
       END INTERFACE

!-----------------------------------------------------------------------
!     SUBROUTINES
!-----------------------------------------------------------------------
       CONTAINS

       LOGICAL FUNCTION adas_tables_avail()
       IMPLICIT NONE
       CHARACTER(LEN=256) :: adasdir, file_str
       adas_tables_avail = .FALSE.
       CALL getenv('ADASDIR', adasdir)
!      Can't INQUIRE about a directory
!       INQUIRE(FILE=TRIM(adasdir),EXIST=adas_tables_avail)
!       IF (.not.adas_tables_avail) THEN
            !WRITE(6,*) 'ADAS FILE: '//TRIM(adasdir)// 'not FOUND!'
!            RETURN
!       END IF
       file_str=TRIM(adasdir) // '/tables/ei/ei_1_coldTarget.cdf'
       INQUIRE(FILE=TRIM(file_str),EXIST=adas_tables_avail)
       IF (.not.adas_tables_avail) THEN
            !WRITE(6,*) 'ADAS FILE: '//TRIM(file_str)// 'not FOUND!'
            RETURN
       END IF
       file_str=TRIM(adasdir) // '/tables/ei/ei_2_coldTarget.cdf'
       INQUIRE(FILE=TRIM(file_str),EXIST=adas_tables_avail)
       IF (.not.adas_tables_avail) THEN
            !WRITE(6,*) 'ADAS FILE: '//TRIM(file_str)// 'not FOUND!'
            RETURN
       END IF
       file_str=TRIM(adasdir) // '/tables/cx/cx_1_1_warmTarget.cdf'
       INQUIRE(FILE=TRIM(file_str),EXIST=adas_tables_avail)
       IF (.not.adas_tables_avail) THEN
            !WRITE(6,*) 'ADAS FILE: '//TRIM(file_str)// 'not FOUND!'
            RETURN
       END IF
       file_str=TRIM(adasdir) // '/tables/cx/cx_2_1_warmTarget.cdf'
       INQUIRE(FILE=TRIM(file_str),EXIST=adas_tables_avail)
       IF (.not.adas_tables_avail) THEN
            !WRITE(6,*) 'ADAS FILE: '//TRIM(file_str)// 'not FOUND!'
            RETURN
       END IF
       file_str=TRIM(adasdir) // '/tables/cx/cx_2_2_warmTarget.cdf'
       INQUIRE(FILE=TRIM(file_str),EXIST=adas_tables_avail)
       IF (.not.adas_tables_avail) THEN
            !WRITE(6,*) 'ADAS FILE: '//TRIM(file_str)// 'not FOUND!'
            RETURN
       END IF
       file_str=TRIM(adasdir) // '/tables/ii/ii_1_1_warmTarget.cdf'
       INQUIRE(FILE=TRIM(file_str),EXIST=adas_tables_avail)
       IF (.not.adas_tables_avail) THEN
            !WRITE(6,*) 'ADAS FILE: '//TRIM(file_str)// 'not FOUND!'
            RETURN
       END IF

       RETURN

       END FUNCTION adas_tables_avail

!       SUBROUTINE adas_load_tables
       SUBROUTINE adas_load_tables(myid,comm)
       IMPLICIT NONE
       INTEGER, INTENT(inout) :: myid
       INTEGER, INTENT(inout) :: comm
       INTEGER, DIMENSION(3) :: dimlen
       CHARACTER(LEN=256) :: adasdir, table_str, var_str

#if defined(MPI_OPT)
       ! Initialize
       CALL deallocate_adas_tables
       ! Load Electron impact iz =1 
       CALL getenv('ADASDIR', adasdir)
       table_str=TRIM(adasdir) // '/tables/ei/ei_1_coldTarget.cdf'
       var_str = 'axis_param_coldTarget'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ei_1_axis, dimlen(1), myid, local_master, comm, win_ei_1a)
       CALL adas_get_var_1D(table_str,var_str,dimlen(1),ei_1_axis,myid,comm)
       var_str = 'sigv'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ei_1_sigv, dimlen(1), myid, local_master, comm, win_ei_1s)
       CALL adas_get_var_1D(table_str,var_str,dimlen(1),ei_1_sigv,myid,comm)

       ! Load Electron impact iz = 2
       table_str=TRIM(adasdir) // '/tables/ei/ei_2_coldTarget.cdf'
       var_str = 'axis_param_coldTarget'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ei_2_axis, dimlen(1), myid, local_master, comm, win_ei_2a)
       CALL adas_get_var_1D(table_str,var_str,dimlen(1),ei_2_axis,myid,comm)
       var_str = 'sigv'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ei_2_sigv, dimlen(1), myid, local_master, comm, win_ei_2s)
       CALL adas_get_var_1D(table_str,var_str,dimlen(1),ei_2_sigv,myid,comm)

       ! Charge Exchange 1-1
       table_str=TRIM(adasdir) // '/tables/cx/cx_1_1_warmTarget.cdf'
       var_str = 'axis_param_warmTarget'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(cx_1_1_axis, dimlen(1), myid, local_master, comm, win_cx_11a)
       CALL adas_get_var_1D(table_str,var_str,dimlen(1),cx_1_1_axis,myid,comm)
       var_str = 'btsigv'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(cx_1_1_btsigv, dimlen(1), dimlen(2), myid, local_master, comm, win_cx_11b)
       CALL adas_get_var_2D(table_str,var_str,dimlen(1),dimlen(2),cx_1_1_btsigv,myid,comm)

       ! Charge Exchange 2-1
       table_str=TRIM(adasdir) // '/tables/cx/cx_2_1_warmTarget.cdf'
       var_str = 'axis_param_warmTarget'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(cx_2_1_axis, dimlen(1), myid, local_master, comm, win_cx_21a)
       CALL adas_get_var_1D(table_str,var_str,dimlen(1),cx_2_1_axis,myid,comm)
       var_str = 'btsigv'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(cx_2_1_btsigv, dimlen(1), dimlen(2), myid, local_master, comm, win_cx_21b)
       CALL adas_get_var_2D(table_str,var_str,dimlen(1),dimlen(2),cx_2_1_btsigv,myid,comm)

       ! Charge Exchange 2-2
       table_str=TRIM(adasdir) // '/tables/cx/cx_2_2_warmTarget.cdf'
       var_str = 'axis_param_warmTarget'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(cx_2_2_axis, dimlen(1), myid, local_master, comm, win_cx_22a)
       CALL adas_get_var_1D(table_str,var_str,dimlen(1),cx_2_2_axis,myid,comm)
       var_str = 'btsigv'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(cx_2_2_btsigv, dimlen(1), dimlen(2), myid, local_master, comm, win_cx_22b)
       CALL adas_get_var_2D(table_str,var_str,dimlen(1),dimlen(2),cx_2_2_btsigv,myid,comm)

       ! Ion Impact 1-1
       table_str=TRIM(adasdir) // '/tables/ii/ii_1_1_warmTarget.cdf'
       var_str = 'axis_param_warmTarget'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ii_1_1_axis, dimlen(1), myid, local_master, comm, win_ii_11a)
       CALL adas_get_var_1D(table_str,var_str,dimlen(1),ii_1_1_axis,myid,comm)
       var_str = 'btsigv'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ii_1_1_btsigv, dimlen(1), dimlen(2), myid, local_master, comm, win_ii_11b)
       CALL adas_get_var_2D(table_str,var_str,dimlen(1),dimlen(2),ii_1_1_btsigv,myid,comm)

       ! Ion Impact 1-2
       table_str=TRIM(adasdir) // '/tables/ii/ii_1_2_warmTarget.cdf'
       var_str = 'axis_param_warmTarget'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ii_1_2_axis, dimlen(1), myid, local_master, comm, win_ii_12a)
       CALL adas_get_var_1D(table_str,var_str,dimlen(1),ii_1_2_axis,myid,comm)
       var_str = 'btsigv'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ii_1_2_btsigv, dimlen(1), dimlen(2), myid, local_master, comm, win_ii_12b)
       CALL adas_get_var_2D(table_str,var_str,dimlen(1),dimlen(2),ii_1_2_btsigv,myid,comm)

       ! Ion Impact 1-3
       table_str=TRIM(adasdir) // '/tables/ii/ii_1_3_warmTarget.cdf'
       var_str = 'axis_param_warmTarget'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ii_1_3_axis, dimlen(1), myid, local_master, comm, win_ii_13a)
       CALL adas_get_var_1D(table_str,var_str,dimlen(1),ii_1_3_axis,myid,comm)
       var_str = 'btsigv'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ii_1_3_btsigv, dimlen(1), dimlen(2), myid, local_master, comm, win_ii_13b)
       CALL adas_get_var_2D(table_str,var_str,dimlen(1),dimlen(2),ii_1_3_btsigv,myid,comm)

       ! Ion Impact 1-4
       table_str=TRIM(adasdir) // '/tables/ii/ii_1_4_warmTarget.cdf'
       var_str = 'axis_param_warmTarget'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ii_1_4_axis, dimlen(1), myid, local_master, comm, win_ii_14a)
       CALL adas_get_var_1D(table_str,var_str,dimlen(1),ii_1_4_axis,myid,comm)
       var_str = 'btsigv'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ii_1_4_btsigv, dimlen(1), dimlen(2), myid, local_master, comm, win_ii_14b)
       CALL adas_get_var_2D(table_str,var_str,dimlen(1),dimlen(2),ii_1_4_btsigv,myid,comm)

       ! Ion Impact 1-5
       table_str=TRIM(adasdir) // '/tables/ii/ii_1_5_warmTarget.cdf'
       var_str = 'axis_param_warmTarget'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ii_1_5_axis, dimlen(1), myid, local_master, comm, win_ii_15a)
       CALL adas_get_var_1D(table_str,var_str,dimlen(1),ii_1_5_axis,myid,comm)
       var_str = 'btsigv'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ii_1_5_btsigv, dimlen(1), dimlen(2), myid, local_master, comm, win_ii_15b)
       CALL adas_get_var_2D(table_str,var_str,dimlen(1),dimlen(2),ii_1_5_btsigv,myid,comm)

       ! Ion Impact 1-6
       table_str=TRIM(adasdir) // '/tables/ii/ii_1_6_warmTarget.cdf'
       var_str = 'axis_param_warmTarget'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ii_1_6_axis, dimlen(1), myid, local_master, comm, win_ii_16a)
       CALL adas_get_var_1D(table_str,var_str,dimlen(1),ii_1_6_axis,myid,comm)
       var_str = 'btsigv'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ii_1_6_btsigv, dimlen(1), dimlen(2), myid, local_master, comm, win_ii_16b)
       CALL adas_get_var_2D(table_str,var_str,dimlen(1),dimlen(2),ii_1_6_btsigv,myid,comm)

       ! Ion Impact 1-7
       table_str=TRIM(adasdir) // '/tables/ii/ii_1_7_warmTarget.cdf'
       var_str = 'axis_param_warmTarget'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ii_1_7_axis, dimlen(1), myid, local_master, comm, win_ii_17a)
       CALL adas_get_var_1D(table_str,var_str,dimlen(1),ii_1_7_axis,myid,comm)
       var_str = 'btsigv'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ii_1_7_btsigv, dimlen(1), dimlen(2), myid, local_master, comm, win_ii_17b)
       CALL adas_get_var_2D(table_str,var_str,dimlen(1),dimlen(2),ii_1_7_btsigv,myid,comm)

       ! Ion Impact 1-8
       table_str=TRIM(adasdir) // '/tables/ii/ii_1_8_warmTarget.cdf'
       var_str = 'axis_param_warmTarget'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ii_1_8_axis, dimlen(1), myid, local_master, comm, win_ii_18a)
       CALL adas_get_var_1D(table_str,var_str,dimlen(1),ii_1_8_axis,myid,comm)
       var_str = 'btsigv'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ii_1_8_btsigv, dimlen(1), dimlen(2), myid, local_master, comm, win_ii_18b)
       CALL adas_get_var_2D(table_str,var_str,dimlen(1),dimlen(2),ii_1_8_btsigv,myid,comm)

       ! Ion Impact 1-9
       table_str=TRIM(adasdir) // '/tables/ii/ii_1_9_warmTarget.cdf'
       var_str = 'axis_param_warmTarget'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ii_1_9_axis, dimlen(1), myid, local_master, comm, win_ii_19a)
       CALL adas_get_var_1D(table_str,var_str,dimlen(1),ii_1_9_axis,myid,comm)
       var_str = 'btsigv'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ii_1_9_btsigv, dimlen(1), dimlen(2), myid, local_master, comm, win_ii_19b)
       CALL adas_get_var_2D(table_str,var_str,dimlen(1),dimlen(2),ii_1_9_btsigv,myid,comm)

       ! Ion Impact 1-10
       table_str=TRIM(adasdir) // '/tables/ii/ii_1_10_warmTarget.cdf'
       var_str = 'axis_param_warmTarget'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ii_1_10_axis, dimlen(1), myid, local_master, comm, win_ii_110a)
       CALL adas_get_var_1D(table_str,var_str,dimlen(1),ii_1_10_axis,myid,comm)
       var_str = 'btsigv'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ii_1_10_btsigv, dimlen(1), dimlen(2), myid, local_master, comm, win_ii_110b)
       CALL adas_get_var_2D(table_str,var_str,dimlen(1),dimlen(2),ii_1_10_btsigv,myid,comm)

       ! Ion Impact 2-1
       table_str=TRIM(adasdir) // '/tables/ii/ii_2_1_warmTarget.cdf'
       var_str = 'axis_param_warmTarget'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ii_2_1_axis, dimlen(1), myid, local_master, comm, win_ii_21a)
       CALL adas_get_var_1D(table_str,var_str,dimlen(1),ii_2_1_axis,myid,comm)
       var_str = 'btsigv'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ii_2_1_btsigv, dimlen(1), dimlen(2), myid, local_master, comm, win_ii_21b)
       CALL adas_get_var_2D(table_str,var_str,dimlen(1),dimlen(2),ii_2_1_btsigv,myid,comm)

       ! Ion Impact 2-2
       table_str=TRIM(adasdir) // '/tables/ii/ii_2_2_warmTarget.cdf'
       var_str = 'axis_param_warmTarget'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ii_2_2_axis, dimlen(1), myid, local_master, comm, win_ii_22a)
       CALL adas_get_var_1D(table_str,var_str,dimlen(1),ii_2_2_axis,myid,comm)
       var_str = 'btsigv'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ii_2_2_btsigv, dimlen(1), dimlen(2), myid, local_master, comm, win_ii_22b)
       CALL adas_get_var_2D(table_str,var_str,dimlen(1),dimlen(2),ii_2_2_btsigv,myid,comm)

       ! Ion Impact 2-4
       table_str=TRIM(adasdir) // '/tables/ii/ii_2_4_warmTarget.cdf'
       var_str = 'axis_param_warmTarget'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ii_2_4_axis, dimlen(1), myid, local_master, comm, win_ii_24a)
       CALL adas_get_var_1D(table_str,var_str,dimlen(1),ii_2_4_axis,myid,comm)
       var_str = 'btsigv'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ii_2_4_btsigv, dimlen(1), dimlen(2), myid, local_master, comm, win_ii_24b)
       CALL adas_get_var_2D(table_str,var_str,dimlen(1),dimlen(2),ii_2_4_btsigv,myid,comm)

       ! Ion Impact 2-6
       table_str=TRIM(adasdir) // '/tables/ii/ii_2_6_warmTarget.cdf'
       var_str = 'axis_param_warmTarget'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ii_2_6_axis, dimlen(1), myid, local_master, comm, win_ii_26a)
       CALL adas_get_var_1D(table_str,var_str,dimlen(1),ii_2_6_axis,myid,comm)
       var_str = 'btsigv'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ii_2_6_btsigv, dimlen(1), dimlen(2), myid, local_master, comm, win_ii_26b)
       CALL adas_get_var_2D(table_str,var_str,dimlen(1),dimlen(2),ii_2_6_btsigv,myid,comm)

       ! Ion Impact 2-8
       table_str=TRIM(adasdir) // '/tables/ii/ii_2_8_warmTarget.cdf'
       var_str = 'axis_param_warmTarget'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ii_2_8_axis, dimlen(1), myid, local_master, comm, win_ii_28a)
       CALL adas_get_var_1D(table_str,var_str,dimlen(1),ii_2_8_axis,myid,comm)
       var_str = 'btsigv'
       CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       CALL mpialloc(ii_2_8_btsigv, dimlen(1), dimlen(2), myid, local_master, comm, win_ii_28b)
       CALL adas_get_var_2D(table_str,var_str,dimlen(1),dimlen(2),ii_2_8_btsigv,myid,comm)

       ! Beam Stopping 1-1
       !table_str=TRIM(adasdir) // '/tables/sv/ex_1_1.cdf'
       !var_str = 'energy_vec'
       !CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       !CALL mpialloc(ex_1_1_e, dimlen(1), myid, local_master, comm, win_ex_11_e)
       !CALL adas_get_var_1D(table_str,var_str,dimlen(1),ex_1_1_e,myid,comm)
       !var_str = 'density_vec'
       !CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       !CALL mpialloc(ex_1_1_n, dimlen(1), myid, local_master, comm, win_ex_11_n)
       !CALL adas_get_var_1D(table_str,var_str,dimlen(1),ex_1_1_n,myid,comm)
       !var_str = 'temp_vec'
       !CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       !CALL mpialloc(ex_1_1_t, dimlen(1), myid, local_master, comm, win_ex_11_t)
       !CALL adas_get_var_1D(table_str,var_str,dimlen(1),ex_1_1_t,myid,comm)
       !var_str = 'sigv_excite'
       !CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       !CALL mpialloc(ex_1_1, dimlen(1), dimlen(2),dimlen(3), myid, local_master, comm, win_ex_11)
       !CALL adas_get_var_3D(table_str,var_str,dimlen(1),dimlen(2),dimlen(3),ex_1_1,myid,comm)

       ! Beam Stopping 1-2
       !table_str=TRIM(adasdir) // '/tables/sv/ex_1_2.cdf'
       !var_str = 'energy_vec'
       !CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       !CALL mpialloc(ex_1_2_e, dimlen(1), myid, local_master, comm, win_ex_12_e)
       !CALL adas_get_var_1D(table_str,var_str,dimlen(1),ex_1_2_e,myid,comm)
       !var_str = 'density_vec'
       !CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       !CALL mpialloc(ex_1_2_n, dimlen(1), myid, local_master, comm, win_ex_12_n)
       !CALL adas_get_var_1D(table_str,var_str,dimlen(1),ex_1_2_n,myid,comm)
       !var_str = 'temp_vec'
       !CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       !CALL mpialloc(ex_1_2_t, dimlen(1), myid, local_master, comm, win_ex_12_t)
       !CALL adas_get_var_1D(table_str,var_str,dimlen(1),ex_1_2_t,myid,comm)
       !var_str = 'sigv_excite'
       !CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       !CALL mpialloc(ex_1_2, dimlen(1), dimlen(2),dimlen(3), myid, local_master, comm, win_ex_12)
       !CALL adas_get_var_3D(table_str,var_str,dimlen(1),dimlen(2),dimlen(3),ex_1_2,myid,comm)

       ! Beam Stopping 1-3
       !table_str=TRIM(adasdir) // '/tables/sv/ex_1_3.cdf'
       !var_str = 'energy_vec'
       !CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       !CALL mpialloc(ex_1_3_e, dimlen(1), myid, local_master, comm, win_ex_13_e)
       !CALL adas_get_var_1D(table_str,var_str,dimlen(1),ex_1_3_e,myid,comm)
       !var_str = 'density_vec'
       !CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       !CALL mpialloc(ex_1_3_n, dimlen(1), myid, local_master, comm, win_ex_13_n)
       !CALL adas_get_var_1D(table_str,var_str,dimlen(1),ex_1_3_n,myid,comm)
       !var_str = 'temp_vec'
       !CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       !CALL mpialloc(ex_1_3_t, dimlen(1), myid, local_master, comm, win_ex_13_t)
       !CALL adas_get_var_1D(table_str,var_str,dimlen(1),ex_1_3_t,myid,comm)
       !var_str = 'sigv_excite'
       !CALL adas_inq_var_size(table_str,var_str,dimlen,myid,comm)
       !CALL mpialloc(ex_1_3, dimlen(1), dimlen(2),dimlen(3), myid, local_master, comm, win_ex_13)
       !CALL adas_get_var_3D(table_str,var_str,dimlen(1),dimlen(2),dimlen(3),ex_1_3,myid,comm)
#endif
       RETURN

       END SUBROUTINE adas_load_tables

       SUBROUTINE deallocate_adas_tables
       IMPLICIT NONE
#if defined(MPI_OPT)
       IF (ASSOCIATED(ei_1_axis)) CALL mpidealloc(ei_1_axis,win_ei_1a)
       IF (ASSOCIATED(ei_2_axis)) CALL mpidealloc(ei_2_axis,win_ei_2a)
       IF (ASSOCIATED(cx_1_1_axis)) CALL mpidealloc(cx_1_1_axis,win_cx_11a)
       IF (ASSOCIATED(cx_2_1_axis)) CALL mpidealloc(cx_2_1_axis,win_cx_21a)
       IF (ASSOCIATED(cx_2_2_axis)) CALL mpidealloc(cx_2_2_axis,win_cx_22a)
       IF (ASSOCIATED(ii_1_1_axis)) CALL mpidealloc(ii_1_1_axis,win_ii_11a)
       IF (ASSOCIATED(ii_1_2_axis)) CALL mpidealloc(ii_1_2_axis,win_ii_12a)
       IF (ASSOCIATED(ii_1_3_axis)) CALL mpidealloc(ii_1_3_axis,win_ii_13a)
       IF (ASSOCIATED(ii_1_4_axis)) CALL mpidealloc(ii_1_4_axis,win_ii_14a)
       IF (ASSOCIATED(ii_1_5_axis)) CALL mpidealloc(ii_1_5_axis,win_ii_15a)
       IF (ASSOCIATED(ii_1_6_axis)) CALL mpidealloc(ii_1_6_axis,win_ii_16a)
       IF (ASSOCIATED(ii_1_7_axis)) CALL mpidealloc(ii_1_7_axis,win_ii_17a)
       IF (ASSOCIATED(ii_1_8_axis)) CALL mpidealloc(ii_1_8_axis,win_ii_18a)
       IF (ASSOCIATED(ii_1_9_axis)) CALL mpidealloc(ii_1_9_axis,win_ii_19a)
       IF (ASSOCIATED(ii_1_10_axis)) CALL mpidealloc(ii_1_10_axis,win_ii_110a)
       IF (ASSOCIATED(ii_2_1_axis)) CALL mpidealloc(ii_2_1_axis,win_ii_21a)
       IF (ASSOCIATED(ii_2_2_axis)) CALL mpidealloc(ii_2_2_axis,win_ii_22a)
       IF (ASSOCIATED(ii_2_3_axis)) CALL mpidealloc(ii_2_3_axis,win_ii_23a)
       IF (ASSOCIATED(ii_2_4_axis)) CALL mpidealloc(ii_2_4_axis,win_ii_24a)
       IF (ASSOCIATED(ii_2_5_axis)) CALL mpidealloc(ii_2_5_axis,win_ii_25a)
       IF (ASSOCIATED(ii_2_6_axis)) CALL mpidealloc(ii_2_6_axis,win_ii_26a)
       IF (ASSOCIATED(ii_2_7_axis)) CALL mpidealloc(ii_2_7_axis,win_ii_27a)
       IF (ASSOCIATED(ii_2_8_axis)) CALL mpidealloc(ii_2_8_axis,win_ii_28a)

       IF (ASSOCIATED(ei_1_sigv)) CALL mpidealloc(ei_1_sigv,win_ei_1s)
       IF (ASSOCIATED(ei_2_sigv)) CALL mpidealloc(ei_2_sigv,win_ei_2s)
       IF (ASSOCIATED(cx_1_1_btsigv)) CALL mpidealloc(cx_1_1_btsigv,win_cx_11b)
       IF (ASSOCIATED(cx_2_1_btsigv)) CALL mpidealloc(cx_2_1_btsigv,win_cx_21b)
       IF (ASSOCIATED(cx_2_2_btsigv)) CALL mpidealloc(cx_2_2_btsigv,win_cx_22b)
       IF (ASSOCIATED(ii_1_1_btsigv)) CALL mpidealloc(ii_1_1_btsigv,win_ii_11b)
       IF (ASSOCIATED(ii_1_2_btsigv)) CALL mpidealloc(ii_1_2_btsigv,win_ii_12b)
       IF (ASSOCIATED(ii_1_3_btsigv)) CALL mpidealloc(ii_1_3_btsigv,win_ii_13b)
       IF (ASSOCIATED(ii_1_4_btsigv)) CALL mpidealloc(ii_1_4_btsigv,win_ii_14b)
       IF (ASSOCIATED(ii_1_5_btsigv)) CALL mpidealloc(ii_1_5_btsigv,win_ii_15b)
       IF (ASSOCIATED(ii_1_6_btsigv)) CALL mpidealloc(ii_1_6_btsigv,win_ii_16b)
       IF (ASSOCIATED(ii_1_7_btsigv)) CALL mpidealloc(ii_1_7_btsigv,win_ii_17b)
       IF (ASSOCIATED(ii_1_8_btsigv)) CALL mpidealloc(ii_1_8_btsigv,win_ii_18b)
       IF (ASSOCIATED(ii_1_9_btsigv)) CALL mpidealloc(ii_1_9_btsigv,win_ii_19b)
       IF (ASSOCIATED(ii_1_10_btsigv)) CALL mpidealloc(ii_1_10_btsigv,win_ii_110b)
       IF (ASSOCIATED(ii_2_1_btsigv)) CALL mpidealloc(ii_2_1_btsigv,win_ii_21b)
       IF (ASSOCIATED(ii_2_2_btsigv)) CALL mpidealloc(ii_2_2_btsigv,win_ii_22b)
       IF (ASSOCIATED(ii_2_3_btsigv)) CALL mpidealloc(ii_2_3_btsigv,win_ii_23b)
       IF (ASSOCIATED(ii_2_4_btsigv)) CALL mpidealloc(ii_2_4_btsigv,win_ii_24b)
       IF (ASSOCIATED(ii_2_5_btsigv)) CALL mpidealloc(ii_2_5_btsigv,win_ii_25b)
       IF (ASSOCIATED(ii_2_6_btsigv)) CALL mpidealloc(ii_2_6_btsigv,win_ii_26b)
       IF (ASSOCIATED(ii_2_7_btsigv)) CALL mpidealloc(ii_2_7_btsigv,win_ii_27b)
       IF (ASSOCIATED(ii_2_8_btsigv)) CALL mpidealloc(ii_2_8_btsigv,win_ii_28b)

       IF (ASSOCIATED(ex_1_1_e)) CALL mpidealloc(ex_1_1_e,win_ex_11_e)
       IF (ASSOCIATED(ex_1_1_n)) CALL mpidealloc(ex_1_1_n,win_ex_11_n)
       IF (ASSOCIATED(ex_1_1_t)) CALL mpidealloc(ex_1_1_t,win_ex_11_t)
       IF (ASSOCIATED(ex_1_2_e)) CALL mpidealloc(ex_1_2_e,win_ex_12_e)
       IF (ASSOCIATED(ex_1_2_n)) CALL mpidealloc(ex_1_2_n,win_ex_12_n)
       IF (ASSOCIATED(ex_1_2_t)) CALL mpidealloc(ex_1_2_t,win_ex_12_t)
       IF (ASSOCIATED(ex_1_3_e)) CALL mpidealloc(ex_1_3_e,win_ex_13_e)
       IF (ASSOCIATED(ex_1_3_n)) CALL mpidealloc(ex_1_3_n,win_ex_13_n)
       IF (ASSOCIATED(ex_1_3_t)) CALL mpidealloc(ex_1_3_t,win_ex_13_t)
       IF (ASSOCIATED(ex_1_1)) CALL mpidealloc(ex_1_1,win_ex_11)
       IF (ASSOCIATED(ex_1_2)) CALL mpidealloc(ex_1_2,win_ex_12)
       IF (ASSOCIATED(ex_1_3)) CALL mpidealloc(ex_1_3,win_ex_13)
#endif
       END SUBROUTINE deallocate_adas_tables

       SUBROUTINE adas_inq_var_size(filename,varname,n,myid,comm)
       IMPLICIT NONE
       CHARACTER(LEN=256), INTENT(in) :: filename, varname
       INTEGER, DIMENSION(3), INTENT(out)  :: n
       INTEGER, INTENT(inout) :: comm, myid
       INTEGER :: ncid, varid, ier, ndims, i
       INTEGER, DIMENSION(3) :: dimid
       n   = 0
#if defined(NETCDF)
       IF (myid == local_master) THEN
          ier = NF90_OPEN(filename,NF90_NOWRITE,ncid)
          ier = NF90_INQ_VARID(ncid,varname,varid)
          ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims, DIMIDS=dimid)
          IF ((ndims) > 3) ndims=3
          DO i = 1, ndims
             ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = n(i))
          END DO
          ier = NF90_CLOSE(ncid)
       END IF
#endif
#if defined(MPI_OPT)
       ier = 0
       CALL MPI_BCAST(n, 2, MPI_INT, local_master, comm, ier)
       RETURN
#endif
       END SUBROUTINE adas_inq_var_size

       SUBROUTINE adas_get_var_1D(filename,varname,n,var,myid,comm)
       IMPLICIT NONE
       CHARACTER(LEN=256), INTENT(in) :: filename, varname
       INTEGER, INTENT(in)  :: n
       DOUBLE PRECISION, DIMENSION(n), INTENT(inout) :: var
       INTEGER, INTENT(inout) :: comm, myid
       INTEGER :: ncid, varid, ier, ndims, i
       INTEGER :: dimid(2)
#if defined(NETCDF)
       IF (myid == local_master) THEN
          ier = NF90_OPEN(filename,NF90_NOWRITE,ncid)
          ier = NF90_INQ_VARID(ncid,varname,varid)
          ier = NF90_GET_VAR(ncid, varid, var)
          ier = NF90_CLOSE(ncid)
       END IF
#endif
#if defined(MPI_OPT)
       ier = 0
       CALL MPI_BCAST(var, n, MPI_DOUBLE, local_master, comm, ier)
       RETURN
#endif
       END SUBROUTINE adas_get_var_1D

       SUBROUTINE adas_get_var_2D(filename,varname,n,m,var,myid,comm)
       IMPLICIT NONE
       CHARACTER(LEN=256), INTENT(in) :: filename, varname
       INTEGER, INTENT(in)  :: n, m
       DOUBLE PRECISION, DIMENSION(n,m), INTENT(inout) :: var
       INTEGER, INTENT(inout) :: comm, myid
       INTEGER :: ncid, varid, ier, ndims, i
       INTEGER :: dimid(2)
#if defined(NETCDF)
       IF (myid == local_master) THEN
          ier = NF90_OPEN(filename,NF90_NOWRITE,ncid)
          ier = NF90_INQ_VARID(ncid,varname,varid)
          ier = NF90_GET_VAR(ncid, varid, var)
          ier = NF90_CLOSE(ncid)
       END IF
#endif
#if defined(MPI_OPT)
       ier = 0
       CALL MPI_BCAST(var, n*m, MPI_DOUBLE, local_master, comm, ier)
       RETURN
#endif
       END SUBROUTINE adas_get_var_2D

       SUBROUTINE adas_get_var_3D(filename,varname,n,m,l,var,myid,comm)
       IMPLICIT NONE
       CHARACTER(LEN=256), INTENT(in) :: filename, varname
       INTEGER, INTENT(in)  :: n, m, l
       DOUBLE PRECISION, DIMENSION(n,m,l), INTENT(inout) :: var
       INTEGER, INTENT(inout) :: comm, myid
       INTEGER :: ncid, varid, ier, ndims, i
       INTEGER :: dimid(3)
#if defined(NETCDF)
       IF (myid == local_master) THEN
          ier = NF90_OPEN(filename,NF90_NOWRITE,ncid)
          ier = NF90_INQ_VARID(ncid,varname,varid)
          ier = NF90_GET_VAR(ncid, varid, var)
          ier = NF90_CLOSE(ncid)
       END IF
#endif
#if defined(MPI_OPT)
       ier = 0
       CALL MPI_BCAST(var, n*m*l, MPI_DOUBLE, local_master, comm, ier)
       RETURN
#endif
       END SUBROUTINE adas_get_var_3D

       SUBROUTINE ADAS_SIGVTE_IONIZ_R8(zneut,tevec,n1,sigv_adas,istat)
       ! calculates <sig*v> electron impact ionization
       real*8, intent(in) :: zneut                   ! atomic charge of primary   particle {1}, available only for H
       real*8, dimension(n1), intent(in) :: tevec      ! vector electron temperature [KeV]
       integer, intent(in) :: n1                       ! vectors size
       real*8, dimension(n1), intent(out) :: sigv_adas !sigma*v [m**3/sec]
       integer, intent(out) :: istat
       integer :: izneut
       izneut=zneut
       CALL adas_sigvte_ioniz_int(izneut,tevec,n1,sigv_adas,istat)
       END SUBROUTINE ADAS_SIGVTE_IONIZ_R8

       SUBROUTINE adas_sigvte_ioniz_int(izneut,tevec,n1,sigv_adas,istat)
       ! calculates <sig*v> electron impact ionization
       integer, intent(in) :: izneut                   ! atomic charge of primary   H and He
       real*8, dimension(n1), intent(in) :: tevec      ! vector electron temperature [KeV]
       integer, intent(in) :: n1                       ! vectors size
       real*8, dimension(n1), intent(out) :: sigv_adas !sigma*v [m**3/sec]
       integer, intent(out) :: istat
       INTEGER :: npts
       REAL*8 :: xlr
       istat = 0
       sigv_adas=0
       IF (izneut == 1) THEN
          xlr = LOG(ei_1_axis(2)/ei_1_axis(1))
          npts = SIZE(ei_1_sigv,DIM=1)
          CALL FLIN1_Z(tevec,sigv_adas,n1,ei_1_sigv,ei_1_axis(1),ei_1_axis(2),xlr,npts,iwarn_adas)
       ELSEIF (izneut == 2) THEN
          xlr = LOG(ei_2_axis(2)/ei_2_axis(1))
          npts = SIZE(ei_2_sigv,DIM=1)
          CALL FLIN1_Z(tevec,sigv_adas,n1,ei_2_sigv,ei_2_axis(1),ei_2_axis(2),xlr,npts,iwarn_adas)
       ELSE
          istat = 1
       END IF
       RETURN
       END SUBROUTINE adas_sigvte_ioniz_int

       SUBROUTINE adas_btsigv_int_int(freact_type,beamchrg,evec,tevec,n1,izneut_in,izion_in,sigv_adas,istat)
       !this subroutine calculates maxwellian average <sig*v> for CX and II reactions
       integer, intent(in) :: freact_type    !reaction type: 1 -- "CX" and 2  -- 'II'
       integer, intent(in) :: izneut_in        !     izneut       : atomic charge of primary   particle {1, 2}
       integer, intent(in) :: izion_in           !     izion        : atomic charge of secondary particle {1, ..., 10}
       INTEGER :: istat  ! completion code; 0 = normal
                                !  istat = 2:  irtype argument invalid
                                !  istat = 3:  izneut or izchrg invalid
                                !  istat = 10:  ADAS data is requested but not found
       real*8, dimension(n1), intent(in) :: evec !vector energy [KeV/amu]
       real*8, dimension(n1), intent(in) :: tevec !vector temperature [KeV]
       integer, intent(in) :: n1  !vectors size
       real*8, dimension(n1), intent(out) :: sigv_adas !sigma*v [m**3/sec]
       integer, intent(in) :: beamchrg ! beam type: =neutral or =ion
       real*8 :: zion_in  
       zion_in=dble(izion_in)
       call adas_btsigv_int_r8(freact_type,beamchrg,evec,tevec,n1,izneut_in,zion_in,sigv_adas,istat)
       END SUBROUTINE adas_btsigv_int_int

       SUBROUTINE adas_btsigv_r8_r8(freact_type,beamchrg,evec,tevec,n1,zneut_in,zion_in,sigv_adas,istat)
       !this subroutine calculates maxwellian average <sig*v>
       integer, intent(in) :: freact_type    !reaction type: 1 -- "CX" and 2  -- 'II'
       real*8, intent(in) :: zneut_in        !     izneut       : atomic charge of primary   particle {1, 2}
       real*8, intent(in) :: zion_in           !     izion        : atomic charge of secondary particle {1, ..., 10}
       INTEGER :: istat  ! completion code; 0 = normal
                                !  istat = 2:  irtype argument invalid
                                !  istat = 3:  izneut or izchrg invalid
                                !  istat = 10:  ADAS data is requested but not found
       real*8, dimension(n1), intent(in) :: evec !vector energy [KeV/amu]
       real*8, dimension(n1), intent(in) :: tevec !vector temperature [KeV]
       integer, intent(in) :: n1  !vectors size
       real*8, dimension(n1), intent(out) :: sigv_adas !sigma*v [m**3/sec]
       integer, intent(in) :: beamchrg ! beam type: =neutral or =ion
       integer :: izneut_in        !     izneut       : atomic charge of primary   particle {1, 2}
       izneut_in=zneut_in
       call adas_btsigv_int_r8(freact_type,beamchrg,evec,tevec,n1,izneut_in,zion_in,sigv_adas,istat)
       END SUBROUTINE adas_btsigv_r8_r8

       SUBROUTINE adas_btsigv_r8_int(freact_type,beamchrg,evec,tevec,n1,zneut_in,izion_in,sigv_adas,istat)
       !this subroutine calculates maxwellian average <sig*v>
       integer, intent(in) :: freact_type    !reaction type: 1 -- "CX" and 2  -- 'II'
       real*8, intent(in) :: zneut_in        !     izneut       : atomic charge of primary   particle {1, 2}
       integer, intent(in) :: izion_in           !     izion        : atomic charge of secondary particle {1, ..., 10}
       INTEGER :: istat  ! completion code; 0 = normal
                                !  istat = 2:  irtype argument invalid
                                !  istat = 3:  izneut or izchrg invalid
                                !  istat = 10:  ADAS data is requested but not found
       real*8, dimension(n1), intent(in) :: evec !vector energy [KeV/amu]
       real*8, dimension(n1), intent(in) :: tevec !vector temperature [KeV]
       integer, intent(in) :: n1  !vectors size
       real*8, dimension(n1), intent(out) :: sigv_adas !sigma*v [m**3/sec]
       integer, intent(in) :: beamchrg ! beam type: =neutral or =ion
       integer :: izneut_in        !     izneut       : atomic charge of primary   particle {1, 2}
       real*8 :: zion_in
       izneut_in=zneut_in
       zion_in=dble(izion_in)
       call adas_btsigv_int_r8(freact_type,beamchrg,evec,tevec,n1,izneut_in,zion_in,sigv_adas,istat)
       END SUBROUTINE adas_btsigv_r8_int

       SUBROUTINE adas_btsigv_int_r8(freact_type,beamchrg,evec,tevec,n1,izneut_in,zion_in,sigv_adas,istat)
       !this subroutine calculates maxwellian average <sig*v>
       integer, intent(in) :: freact_type    !reaction type: 1 -- "CX" and 2  -- 'II'
       integer, intent(in) :: izneut_in        !     izneut       : atomic charge of primary   particle {1, 2}
       real*8, intent(in) :: zion_in           !     izion        : atomic charge of secondary particle {1, ..., 10}
       INTEGER :: istat  ! completion code; 0 = normal
                                !  istat = 2:  irtype argument invalid
                                !  istat = 3:  izneut or izchrg invalid
                                !  istat = 10:  ADAS data is requested but not found
       real*8, dimension(n1), intent(in) :: evec !vector energy [KeV/amu]
       real*8, dimension(n1), intent(in) :: tevec !vector temperature [KeV]
       integer, intent(in) :: n1  !vectors size
       real*8, dimension(n1), intent(out) :: sigv_adas !sigma*v [m**3/sec]
       integer, intent(in) :: beamchrg ! beam type: =neutral or =ion

       integer, parameter :: NEUTRAL=1, ION=2
       real*8  :: zbeam         !     zbeam        : atomic charge of primate particle
       real*8  :: zion          !     izion        : atomic charge of secondary particle {1, ..., 10}
       integer :: npts_e, npts_t ! number of points for table
       real*8  :: xlr, xlr_t ! If XLR>0, then the X grid is equally spaced on a logarithmic scale:
                  ! LINEAR IF XLR .LE. 0.0
       real*8, dimension(n1) :: sigv_adas_wrk1,sigv_adas_wrk2 !sigma*v [m**3/sec]
       !checking input data for availability
       istat=0
       sigv_adas=0
       sigv_adas_wrk1 =0; sigv_adas_wrk2=0

       ! The the beam is an ion then flip it.
       IF (beamchrg.eq.NEUTRAL) THEN
          zbeam = DBLE(izneut_in)
          zion  = DBLE(zion_in)
       ELSE ! FLIP for ION-BEAM
          zbeam = DBLE(zion_in)
          zion  = DBLE(izneut_in)
       END IF

       IF (freact_type == 1) THEN ! Charge Exchange (cx_zneut_zion)
          IF (zbeam == 1) THEN
             xlr    =  log(cx_1_1_axis(2)/cx_1_1_axis(1))
             xlr_t  =  log(cx_1_1_axis(6)/cx_1_1_axis(5))
             npts_e = SIZE(cx_1_1_btsigv,DIM=1)
             npts_t = SIZE(cx_1_1_btsigv,DIM=2)
             call FLINT_Z(evec,tevec,sigv_adas_wrk1,n1,cx_1_1_btsigv, cx_1_1_axis(1),cx_1_1_axis(2),&
                                               xlr,npts_e,cx_1_1_axis(5),cx_1_1_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
             sigv_adas = sigv_adas_wrk1*(zion) ! This is the ADAS WAY
          ELSE IF (zbeam == 2) THEN
             IF (zion >=1 .and. zion < 2) THEN
                xlr=log(cx_2_1_axis(2)/cx_2_1_axis(1))
                xlr_t=log(cx_2_1_axis(6)/cx_2_1_axis(5))
                npts_e = SIZE(cx_2_1_btsigv,DIM=1)
                npts_t = SIZE(cx_2_1_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk1,n1,cx_2_1_btsigv,cx_2_1_axis(1),cx_2_1_axis(2),&
                             xlr,npts_e,cx_2_1_axis(5),cx_2_1_axis(6),&
                             xlr_t,npts_t,iwarn_adas)
                xlr=log(cx_2_2_axis(2)/cx_2_2_axis(1))
                xlr_t=log(cx_2_2_axis(6)/cx_2_2_axis(5))
                npts_e = SIZE(cx_2_2_btsigv,DIM=1)
                npts_t = SIZE(cx_2_2_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk2,n1,cx_2_2_btsigv,cx_2_2_axis(1),cx_2_2_axis(2),&
                             xlr,npts_e,cx_2_2_axis(5),cx_2_2_axis(6),&
                             xlr_t,npts_t,iwarn_adas)
                sigv_adas = sigv_adas_wrk1 + (sigv_adas_wrk2-sigv_adas_wrk1)*(zion-1)
             ELSE IF (zion >=2) THEN
                xlr    =  log(cx_2_2_axis(2)/cx_2_2_axis(1))
                xlr_t  =  log(cx_2_2_axis(6)/cx_2_2_axis(5))
                npts_e = SIZE(cx_2_2_btsigv,DIM=1)
                npts_t = SIZE(cx_2_2_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk1,n1,cx_2_2_btsigv, cx_2_2_axis(1),cx_2_2_axis(2),&
                                               xlr,npts_e,cx_2_2_axis(5),cx_2_2_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
                sigv_adas = sigv_adas_wrk1*(zion/10) ! This is the ADAS WAY
             END IF
          ELSE
             istat = 1
          END IF
       ELSE IF (freact_type == 2) THEN
          IF(zbeam == 1) THEN
             IF (zion >=1 .and. zion < 2) THEN
                xlr=log(ii_1_1_axis(2)/ii_1_1_axis(1))
                xlr_t=log(ii_1_1_axis(6)/ii_1_1_axis(5))
                npts_e = SIZE(ii_1_1_btsigv,DIM=1)
                npts_t = SIZE(ii_1_1_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk1,n1,ii_1_1_btsigv,ii_1_1_axis(1),ii_1_1_axis(2),&
                             xlr,npts_e,ii_1_1_axis(5),ii_1_1_axis(6),&
                             xlr_t,npts_t,iwarn_adas)
                xlr=log(ii_1_2_axis(2)/ii_1_2_axis(1))
                xlr_t=log(ii_1_2_axis(6)/ii_1_2_axis(5))
                npts_e = SIZE(ii_1_2_btsigv,DIM=1)
                npts_t = SIZE(ii_1_2_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk2,n1,ii_1_2_btsigv,ii_1_2_axis(1),ii_1_2_axis(2),&
                             xlr,npts_e,ii_1_2_axis(5),ii_1_2_axis(6),&
                             xlr_t,npts_t,iwarn_adas)
                sigv_adas = sigv_adas_wrk1 + (sigv_adas_wrk2-sigv_adas_wrk1)*(zion-1)
             ELSE IF (zion >=2 .and. zion < 3) THEN
                xlr=log(ii_1_2_axis(2)/ii_1_2_axis(1))
                xlr_t=log(ii_1_2_axis(6)/ii_1_2_axis(5))
                npts_e = SIZE(ii_1_2_btsigv,DIM=1)
                npts_t = SIZE(ii_1_2_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk1,n1,ii_1_2_btsigv,ii_1_2_axis(1),ii_1_2_axis(2),&
                             xlr,npts_e,ii_1_2_axis(5),ii_1_2_axis(6),&
                             xlr_t,npts_t,iwarn_adas)
                xlr=log(ii_1_3_axis(2)/ii_1_3_axis(1))
                xlr_t=log(ii_1_3_axis(6)/ii_1_3_axis(5))
                npts_e = SIZE(ii_1_3_btsigv,DIM=1)
                npts_t = SIZE(ii_1_3_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk2,n1,ii_1_3_btsigv,ii_1_3_axis(1),ii_1_3_axis(2),&
                             xlr,npts_e,ii_1_3_axis(5),ii_1_3_axis(6),&
                             xlr_t,npts_t,iwarn_adas)
                sigv_adas = sigv_adas_wrk1 + (sigv_adas_wrk2-sigv_adas_wrk1)*(zion-2)
             ELSE IF (zion >=3 .and. zion < 4) THEN
                xlr    =  log(ii_1_3_axis(2)/ii_1_3_axis(1))
                xlr_t  =  log(ii_1_3_axis(6)/ii_1_3_axis(5))
                npts_e = SIZE(ii_1_3_btsigv,DIM=1)
                npts_t = SIZE(ii_1_3_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk1,n1,ii_1_3_btsigv, ii_1_3_axis(1),ii_1_3_axis(2),&
                                               xlr,npts_e,ii_1_3_axis(5),ii_1_3_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
                xlr    =  log(ii_1_4_axis(2)/ii_1_4_axis(1))
                xlr_t  =  log(ii_1_4_axis(6)/ii_1_4_axis(5))
                npts_e = SIZE(ii_1_4_btsigv,DIM=1)
                npts_t = SIZE(ii_1_4_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk2,n1,ii_1_4_btsigv, ii_1_4_axis(1),ii_1_4_axis(2),&
                                               xlr,npts_e,ii_1_4_axis(5),ii_1_4_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
                sigv_adas = sigv_adas_wrk1 + (sigv_adas_wrk2-sigv_adas_wrk1)*(zion-3)
             ELSE IF (zion >=4 .and. zion < 5) THEN
                xlr    =  log(ii_1_4_axis(2)/ii_1_4_axis(1))
                xlr_t  =  log(ii_1_4_axis(6)/ii_1_4_axis(5))
                npts_e = SIZE(ii_1_4_btsigv,DIM=1)
                npts_t = SIZE(ii_1_4_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk1,n1,ii_1_4_btsigv, ii_1_4_axis(1),ii_1_4_axis(2),&
                                               xlr,npts_e,ii_1_4_axis(5),ii_1_4_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
                xlr    =  log(ii_1_5_axis(2)/ii_1_5_axis(1))
                xlr_t  =  log(ii_1_5_axis(6)/ii_1_5_axis(5))
                npts_e = SIZE(ii_1_5_btsigv,DIM=1)
                npts_t = SIZE(ii_1_5_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk2,n1,ii_1_5_btsigv, ii_1_5_axis(1),ii_1_5_axis(2),&
                                               xlr,npts_e,ii_1_5_axis(5),ii_1_5_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
                sigv_adas = sigv_adas_wrk1 + (sigv_adas_wrk2-sigv_adas_wrk1)*(zion-4)
             ELSE IF (zion >=5 .and. zion < 6) THEN
                xlr    =  log(ii_1_5_axis(2)/ii_1_5_axis(1))
                xlr_t  =  log(ii_1_5_axis(6)/ii_1_5_axis(5))
                npts_e = SIZE(ii_1_5_btsigv,DIM=1)
                npts_t = SIZE(ii_1_5_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk1,n1,ii_1_5_btsigv, ii_1_5_axis(1),ii_1_5_axis(2),&
                                               xlr,npts_e,ii_1_5_axis(5),ii_1_5_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
                xlr    =  log(ii_1_6_axis(2)/ii_1_6_axis(1))
                xlr_t  =  log(ii_1_6_axis(6)/ii_1_6_axis(5))
                npts_e = SIZE(ii_1_6_btsigv,DIM=1)
                npts_t = SIZE(ii_1_6_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk2,n1,ii_1_6_btsigv, ii_1_6_axis(1),ii_1_6_axis(2),&
                                               xlr,npts_e,ii_1_6_axis(5),ii_1_6_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
                sigv_adas = sigv_adas_wrk1 + (sigv_adas_wrk2-sigv_adas_wrk1)*(zion-5)
             ELSE IF (zion >=6 .and. zion < 7) THEN
                xlr    =  log(ii_1_6_axis(2)/ii_1_6_axis(1))
                xlr_t  =  log(ii_1_6_axis(6)/ii_1_6_axis(5))
                npts_e = SIZE(ii_1_6_btsigv,DIM=1)
                npts_t = SIZE(ii_1_6_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk1,n1,ii_1_6_btsigv, ii_1_6_axis(1),ii_1_6_axis(2),&
                                               xlr,npts_e,ii_1_6_axis(5),ii_1_6_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
                xlr    =  log(ii_1_7_axis(2)/ii_1_7_axis(1))
                xlr_t  =  log(ii_1_7_axis(6)/ii_1_7_axis(5))
                npts_e = SIZE(ii_1_7_btsigv,DIM=1)
                npts_t = SIZE(ii_1_7_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk2,n1,ii_1_7_btsigv, ii_1_7_axis(1),ii_1_7_axis(2),&
                                               xlr,npts_e,ii_1_7_axis(5),ii_1_7_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
                sigv_adas = sigv_adas_wrk1 + (sigv_adas_wrk2-sigv_adas_wrk1)*(zion-6)
             ELSE IF (zion >=7 .and. zion < 8) THEN
                xlr    =  log(ii_1_7_axis(2)/ii_1_7_axis(1))
                xlr_t  =  log(ii_1_7_axis(6)/ii_1_7_axis(5))
                npts_e = SIZE(ii_1_7_btsigv,DIM=1)
                npts_t = SIZE(ii_1_7_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk1,n1,ii_1_7_btsigv, ii_1_7_axis(1),ii_1_7_axis(2),&
                                               xlr,npts_e,ii_1_7_axis(5),ii_1_7_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
                xlr    =  log(ii_1_8_axis(2)/ii_1_8_axis(1))
                xlr_t  =  log(ii_1_8_axis(6)/ii_1_8_axis(5))
                npts_e = SIZE(ii_1_8_btsigv,DIM=1)
                npts_t = SIZE(ii_1_8_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk2,n1,ii_1_8_btsigv, ii_1_8_axis(1),ii_1_8_axis(2),&
                                               xlr,npts_e,ii_1_8_axis(5),ii_1_8_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
                sigv_adas = sigv_adas_wrk1 + (sigv_adas_wrk2-sigv_adas_wrk1)*(zion-7)
             ELSE IF (zion >=8 .and. zion < 9) THEN
                xlr    =  log(ii_1_8_axis(2)/ii_1_8_axis(1))
                xlr_t  =  log(ii_1_8_axis(6)/ii_1_8_axis(5))
                npts_e = SIZE(ii_1_8_btsigv,DIM=1)
                npts_t = SIZE(ii_1_8_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk1,n1,ii_1_8_btsigv, ii_1_8_axis(1),ii_1_8_axis(2),&
                                               xlr,npts_e,ii_1_8_axis(5),ii_1_8_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
                xlr    =  log(ii_1_9_axis(2)/ii_1_9_axis(1))
                xlr_t  =  log(ii_1_9_axis(6)/ii_1_9_axis(5))
                npts_e = SIZE(ii_1_9_btsigv,DIM=1)
                npts_t = SIZE(ii_1_9_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk2,n1,ii_1_9_btsigv, ii_1_9_axis(1),ii_1_9_axis(2),&
                                               xlr,npts_e,ii_1_9_axis(5),ii_1_9_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
                sigv_adas = sigv_adas_wrk1 + (sigv_adas_wrk2-sigv_adas_wrk1)*(zion-8)
             ELSE IF (zion >=9 .and. zion < 10) THEN
                xlr    =  log(ii_1_9_axis(2)/ii_1_9_axis(1))
                xlr_t  =  log(ii_1_9_axis(6)/ii_1_9_axis(5))
                npts_e = SIZE(ii_1_9_btsigv,DIM=1)
                npts_t = SIZE(ii_1_9_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk1,n1,ii_1_9_btsigv, ii_1_9_axis(1),ii_1_9_axis(2),&
                                               xlr,npts_e,ii_1_9_axis(5),ii_1_9_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
                xlr    =  log(ii_1_10_axis(2)/ii_1_10_axis(1))
                xlr_t  =  log(ii_1_10_axis(6)/ii_1_10_axis(5))
                npts_e = SIZE(ii_1_10_btsigv,DIM=1)
                npts_t = SIZE(ii_1_10_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk2,n1,ii_1_10_btsigv, ii_1_10_axis(1),ii_1_10_axis(2),&
                                               xlr,npts_e,ii_1_10_axis(5),ii_1_10_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
                sigv_adas = sigv_adas_wrk1 + (sigv_adas_wrk2-sigv_adas_wrk1)*(zion-9)
             ELSE IF (zion >=10) THEN
                xlr    =  log(ii_1_10_axis(2)/ii_1_10_axis(1))
                xlr_t  =  log(ii_1_10_axis(6)/ii_1_10_axis(5))
                npts_e = SIZE(ii_1_10_btsigv,DIM=1)
                npts_t = SIZE(ii_1_10_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk1,n1,ii_1_10_btsigv, ii_1_10_axis(1),ii_1_10_axis(2),&
                                               xlr,npts_e,ii_1_10_axis(5),ii_1_10_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
                sigv_adas = sigv_adas_wrk1*(zion/10) ! This is the ADAS WAY
             END IF
          ELSE IF (zbeam == 2) THEN
             IF (zion >=1 .and. zion < 2) THEN
                xlr=log(ii_2_1_axis(2)/ii_2_1_axis(1))
                xlr_t=log(ii_2_1_axis(6)/ii_2_1_axis(5))
                npts_e = SIZE(ii_2_1_btsigv,DIM=1)
                npts_t = SIZE(ii_2_1_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk1,n1,ii_2_1_btsigv,ii_2_1_axis(1),ii_2_1_axis(2),&
                             xlr,npts_e,ii_2_1_axis(5),ii_2_1_axis(6),&
                             xlr_t,npts_t,iwarn_adas)
                xlr=log(ii_2_2_axis(2)/ii_2_2_axis(1))
                xlr_t=log(ii_2_2_axis(6)/ii_2_2_axis(5))
                npts_e = SIZE(ii_2_2_btsigv,DIM=1)
                npts_t = SIZE(ii_2_2_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk2,n1,ii_2_2_btsigv,ii_2_2_axis(1),ii_2_2_axis(2),&
                             xlr,npts_e,ii_2_2_axis(5),ii_2_2_axis(6),&
                             xlr_t,npts_t,iwarn_adas)
                sigv_adas = sigv_adas_wrk1 + (sigv_adas_wrk2-sigv_adas_wrk1)*(zion-1)
             ELSE IF (zion >=2 .and. zion < 4) THEN
                xlr=log(ii_2_2_axis(2)/ii_2_2_axis(1))
                xlr_t=log(ii_2_2_axis(6)/ii_2_2_axis(5))
                npts_e = SIZE(ii_2_2_btsigv,DIM=1)
                npts_t = SIZE(ii_2_2_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk1,n1,ii_2_2_btsigv,ii_2_2_axis(1),ii_2_2_axis(2),&
                             xlr,npts_e,ii_2_2_axis(5),ii_2_2_axis(6),&
                             xlr_t,npts_t,iwarn_adas)
                xlr    =  log(ii_2_4_axis(2)/ii_2_4_axis(1))
                xlr_t  =  log(ii_2_4_axis(6)/ii_2_4_axis(5))
                npts_e = SIZE(ii_2_4_btsigv,DIM=1)
                npts_t = SIZE(ii_2_4_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk2,n1,ii_2_4_btsigv,ii_2_4_axis(1),ii_2_4_axis(2),&
                             xlr,npts_e,ii_2_4_axis(5),ii_2_4_axis(6),&
                             xlr_t,npts_t,iwarn_adas)
                sigv_adas = sigv_adas_wrk1 + (sigv_adas_wrk2-sigv_adas_wrk1)*(zion-2)
             ELSE IF (zion >=4 .and. zion < 6) THEN
                xlr    =  log(ii_2_4_axis(2)/ii_2_4_axis(1))
                xlr_t  =  log(ii_2_4_axis(6)/ii_2_4_axis(5))
                npts_e = SIZE(ii_2_4_btsigv,DIM=1)
                npts_t = SIZE(ii_2_4_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk1,n1,ii_2_4_btsigv, ii_2_4_axis(1),ii_2_4_axis(2),&
                                               xlr,npts_e,ii_2_4_axis(5),ii_2_4_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
                xlr    =  log(ii_2_6_axis(2)/ii_2_6_axis(1))
                xlr_t  =  log(ii_2_6_axis(6)/ii_2_6_axis(5))
                npts_e = SIZE(ii_2_6_btsigv,DIM=1)
                npts_t = SIZE(ii_2_6_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk2,n1,ii_2_6_btsigv, ii_2_6_axis(1),ii_2_6_axis(2),&
                                               xlr,npts_e,ii_2_6_axis(5),ii_2_6_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
                sigv_adas = sigv_adas_wrk1 + (sigv_adas_wrk2-sigv_adas_wrk1)*(zion-4)
             ELSE IF (zion >=6 .and. zion < 8) THEN
                xlr    =  log(ii_2_6_axis(2)/ii_2_6_axis(1))
                xlr_t  =  log(ii_2_6_axis(6)/ii_2_6_axis(5))
                npts_e = SIZE(ii_2_6_btsigv,DIM=1)
                npts_t = SIZE(ii_2_6_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk1,n1,ii_2_6_btsigv, ii_2_6_axis(1),ii_2_6_axis(2),&
                                               xlr,npts_e,ii_2_6_axis(5),ii_2_6_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
                xlr    =  log(ii_2_8_axis(2)/ii_2_8_axis(1))
                xlr_t  =  log(ii_2_8_axis(6)/ii_2_8_axis(5))
                npts_e = SIZE(ii_2_8_btsigv,DIM=1)
                npts_t = SIZE(ii_2_8_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk2,n1,ii_2_8_btsigv, ii_2_8_axis(1),ii_2_8_axis(2),&
                                               xlr,npts_e,ii_2_8_axis(5),ii_2_8_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
                sigv_adas = sigv_adas_wrk1 + (sigv_adas_wrk2-sigv_adas_wrk1)*(zion-6)
             ELSE IF (zion >=8) THEN
                xlr    =  log(ii_2_8_axis(2)/ii_2_8_axis(1))
                xlr_t  =  log(ii_2_8_axis(6)/ii_2_8_axis(5))
                npts_e = SIZE(ii_2_8_btsigv,DIM=1)
                npts_t = SIZE(ii_2_8_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk1,n1,ii_2_8_btsigv, ii_2_8_axis(1),ii_2_8_axis(2),&
                                               xlr,npts_e,ii_2_8_axis(5),ii_2_8_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
                sigv_adas = sigv_adas_wrk1*(zion/8) ! This is the ADAS WAY
             END IF
          ELSE
             istat = 1
          END IF
       END IF
       RETURN
       END SUBROUTINE adas_btsigv_int_r8

       SUBROUTINE adas_bms_int_r8(EB, Nevec, ZB, CONC, TI, ZI, NE, SIGVI, iERR)
       !Compute beam stopping rate coefficient on impurity
!                            9 = array initialisation error)
!
       REAL*8,  INTENT(IN)  :: EB(nevec)    ! Beam energy [KeV/amu]
       INTEGER, INTENT(IN)  :: nevec        ! Energy vector size
       INTEGER, INTENT(IN)  :: ZB           ! Beam Particle Charge Number
       REAL*8,  INTENT(IN)  :: conc(nevec)  !target species fraction
!                                            : Z_imp * N_imp/SUM(Z_imp*N_imp) (0 < conc(i) <= 1)
!                                            : where,
!                                            : Z_imp -- Particle Nuclear Charge Number 
!                                            : N_imp -- impurity density
       REAL*8,  INTENT(IN)  :: TI(nevec)    ! Target temperature [KeV/amu] 
       REAL*8,  INTENT(IN)  :: ZI           ! Target Nuclear Charge Number
       REAL*8,  INTENT(IN)  :: NE(nevec)    ! Equivalent electron density [cm**-3]
!                                              N_el*SUM(Z_imp^2* N_imp)/SUM(Z_imp * N_imp)/Z_imp
       REAL*8,  INTENT(OUT) :: SIGVI(nevec) !Beam Stopping Rate Coefficient [m**3/s]
       INTEGER              :: iERR
!     iERR     : Error Code (0 = OK,                             
!                            1 = unsupported beam species
!                            2 = unsupported target species
!                            3 = input data corrupted (negative beam energy, densities or temperatures)
!                            4 = input data out of range
!                            5 = adas file reading error
!                            6 = spline output data invalid (negative)
!                            7 = spline error
!                            8 = array allocation error
!                            9 = array initialisation error)
       INTEGER :: n, m, l, i
       REAL*8  :: nfrac, mfrac, lfrac, sigvi_wrk

       SIGVI = 0
       iERR  = 0
       IF (ZB == 1) THEN ! H beam
            IF (ZI < 1) THEN
               PRINT *,'ERROR: ZI_r8 < 1'
               iERR = 2
               RETURN
            ELSE IF (ZI <2) THEN
               DO i=1, nevec
                  n = MIN(COUNT(ex_1_1_e < EB(i)),1)
                  m = MIN(COUNT(ex_1_1_n < NE(i)),1)
                  l = MIN(COUNT(ex_1_1_t < TI(i)),1)
                  nfrac = (EB(i) - ex_1_1_e(n))/(ex_1_1_e(n+1)-ex_1_1_e(n))
                  mfrac = (NE(i) - ex_1_1_n(m))/(ex_1_1_n(m+1)-ex_1_1_n(m))
                  lfrac = (TI(i) - ex_1_1_t(l))/(ex_1_1_t(l+1)-ex_1_1_t(l))
                  sigvi(i) = ex_1_1(n,m,l) + nfrac*(ex_1_1(n+1,m  ,l  )-ex_1_1(n,m,l)) &
                                            + mfrac*(ex_1_1(n  ,m+1,l  )-ex_1_1(n,m,l)) &
                                            + lfrac*(ex_1_1(n  ,m  ,l+1)-ex_1_1(n,m,l))
                  n = MIN(COUNT(ex_1_2_e < EB(i)),1)
                  m = MIN(COUNT(ex_1_2_n < NE(i)),1)
                  l = MIN(COUNT(ex_1_2_t < TI(i)),1)
                  nfrac = (EB(i) - ex_1_2_e(n))/(ex_1_2_e(n+1)-ex_1_2_e(n))
                  mfrac = (NE(i) - ex_1_2_n(m))/(ex_1_2_n(m+1)-ex_1_2_n(m))
                  lfrac = (TI(i) - ex_1_2_t(l))/(ex_1_2_t(l+1)-ex_1_2_t(l))
                  sigvi_wrk = ex_1_2(n,m,l) + nfrac*(ex_1_2(n+1,m  ,l  )-ex_1_2(n,m,l)) &
                                            + mfrac*(ex_1_2(n  ,m+1,l  )-ex_1_2(n,m,l)) &
                                            + lfrac*(ex_1_2(n  ,m  ,l+1)-ex_1_2(n,m,l))
                  sigvi(i) = sigvi(i) + (sigvi_wrk-sigvi(i)) * (ZI-1)
               END DO
            ELSE IF (ZI <3) THEN
               DO i=1, nevec
                  n = MIN(COUNT(ex_1_2_e < EB(i)),1)
                  m = MIN(COUNT(ex_1_2_n < NE(i)),1)
                  l = MIN(COUNT(ex_1_2_t < TI(i)),1)
                  nfrac = (EB(i) - ex_1_2_e(n))/(ex_1_2_e(n+1)-ex_1_2_e(n))
                  mfrac = (NE(i) - ex_1_2_n(m))/(ex_1_2_n(m+1)-ex_1_2_n(m))
                  lfrac = (TI(i) - ex_1_2_t(l))/(ex_1_2_t(l+1)-ex_1_2_t(l))
                  sigvi(i) = ex_1_2(n,m,l) + nfrac*(ex_1_2(n+1,m  ,l  )-ex_1_2(n,m,l)) &
                                            + mfrac*(ex_1_2(n  ,m+1,l  )-ex_1_2(n,m,l)) &
                                            + lfrac*(ex_1_2(n  ,m  ,l+1)-ex_1_2(n,m,l))
                  n = MIN(COUNT(ex_1_3_e < EB(i)),1)
                  m = MIN(COUNT(ex_1_3_n < NE(i)),1)
                  l = MIN(COUNT(ex_1_3_t < TI(i)),1)
                  nfrac = (EB(i) - ex_1_3_e(n))/(ex_1_3_e(n+1)-ex_1_3_e(n))
                  mfrac = (NE(i) - ex_1_3_n(m))/(ex_1_3_n(m+1)-ex_1_3_n(m))
                  lfrac = (TI(i) - ex_1_3_t(l))/(ex_1_3_t(l+1)-ex_1_3_t(l))
                  sigvi_wrk = ex_1_3(n,m,l) + nfrac*(ex_1_3(n+1,m  ,l  )-ex_1_3(n,m,l)) &
                                            + mfrac*(ex_1_3(n  ,m+1,l  )-ex_1_3(n,m,l)) &
                                            + lfrac*(ex_1_3(n  ,m  ,l+1)-ex_1_3(n,m,l))
                  sigvi(i) = sigvi(i) + (sigvi_wrk-sigvi(i)) * (ZI-2)
               END DO
            END IF

       ELSE
            PRINT *,'ERROR: Non-Hydrogenic species missing'
            iERR = 1
       END IF
       RETURN
       END SUBROUTINE adas_bms_int_r8

end module adas_mod_parallel