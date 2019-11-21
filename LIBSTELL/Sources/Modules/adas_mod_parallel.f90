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
!       USE mpi_sharmem
#endif
!-----------------------------------------------------------------------
!     Variables
!-----------------------------------------------------------------------
       IMPLICIT NONE
       INTEGER :: iwarn_adas
       DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: ei_1_axis, ei_2_axis
       DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: ei_1_sigv, ei_2_sigv
       DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: cx_1_1_axis, cx_2_1_axis,cx_2_2_axis,&
                                                      ii_1_1_axis, ii_1_2_axis, &
                                                      ii_1_3_axis, ii_1_4_axis, &
                                                      ii_1_5_axis, ii_1_6_axis, &
                                                      ii_1_7_axis, ii_1_8_axis, &
                                                      ii_1_9_axis, ii_1_10_axis, &
                                                      ii_2_1_axis, ii_2_2_axis, &
                                                      ii_2_3_axis, ii_2_4_axis, &
                                                      ii_2_5_axis, ii_2_6_axis, &
                                                      ii_2_7_axis, ii_2_8_axis
       DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: cx_1_1_btsigv,cx_2_1_btsigv,cx_2_2_btsigv,&
                                                        ii_1_1_btsigv, ii_1_2_btsigv, &
                                                        ii_1_3_btsigv, ii_1_4_btsigv, &
                                                        ii_1_5_btsigv, ii_1_6_btsigv, &
                                                        ii_1_7_btsigv, ii_1_8_btsigv, &
                                                        ii_1_9_btsigv, ii_1_10_btsigv, &
                                                        ii_2_1_btsigv, ii_2_2_btsigv, &
                                                        ii_2_3_btsigv, ii_2_4_btsigv, &
                                                        ii_2_5_btsigv, ii_2_6_btsigv, &
                                                        ii_2_7_btsigv, ii_2_8_btsigv
!-----------------------------------------------------------------------
!     INTERFACES
!-----------------------------------------------------------------------
       PUBLIC :: adas_btsigv
       PUBLIC :: adas_sigvte_ioniz
       
       INTERFACE adas_sigvte_ioniz
          MODULE PROCEDURE adas_sigvte_ioniz_r8,adas_sigvte_ioniz_int
       END INTERFACE

       INTERFACE adas_btsigv
          MODULE PROCEDURE adas_btsigv_int_int,adas_btsigv_int_r8,adas_btsigv_r8_r8,adas_btsigv_r8_int
       END INTERFACE

!-----------------------------------------------------------------------
!     SUBROUTINES
!-----------------------------------------------------------------------
       CONTAINS

       SUBROUTINE adas_load_tables
       IMPLICIT NONE
       INTEGER :: ncid, varid, ier, ndims, i
       INTEGER, ALLOCATABLE, DIMENSION(:) :: dimid, dimlen
       CHARACTER(LEN=256) :: adasdir, table_str

       ! Initialize
       CALL deallocate_adas_tables

#if defined(NETCDF)
       ! Load Electron impact iz =1 
       CALL getenv('ADASDIR', adasdir)
       table_str=TRIM(adasdir) // '/tables/ei/ei_1_coldTarget.cdf'
       ier = NF90_OPEN(table_str,NF90_NOWRITE,ncid)
       ier = NF90_INQ_VARID(ncid,'axis_param_coldTarget',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ei_1_axis(dimlen(1)))
       ier = NF90_GET_VAR(ncid, varid, ei_1_axis)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_INQ_VARID(ncid,'sigv',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ei_1_sigv(dimlen(1)))
       ier = NF90_GET_VAR(ncid, varid, ei_1_sigv)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_CLOSE(ncid)

       ! Load Electron impact iz =2 
       table_str=TRIM(adasdir) // '/tables/ei/ei_2_coldTarget.cdf'
       ier = NF90_OPEN(table_str,NF90_NOWRITE,ncid)
       ier = NF90_INQ_VARID(ncid,'axis_param_coldTarget',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ei_2_axis(dimlen(1)))
       ier = NF90_GET_VAR(ncid, varid, ei_2_axis)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_INQ_VARID(ncid,'sigv',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ei_2_sigv(dimlen(1)))
       ier = NF90_GET_VAR(ncid, varid, ei_2_sigv)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_CLOSE(ncid)

       ! Load Ion Impact CX 1-1
       table_str=TRIM(adasdir) // '/tables/cx/cx_1_1_warmTarget.cdf'
       ier = NF90_OPEN(table_str,NF90_NOWRITE,ncid)
       ier = NF90_INQ_VARID(ncid,'axis_param_warmTarget',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(cx_1_1_axis(dimlen(1)))
       ier = NF90_GET_VAR(ncid, varid, cx_1_1_axis)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_INQ_VARID(ncid,'btsigv',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(cx_1_1_btsigv(dimlen(1),dimlen(2)))
       ier = NF90_GET_VAR(ncid, varid, cx_1_1_btsigv)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_CLOSE(ncid)

       ! Load Ion Impact CX 2-1
       table_str=TRIM(adasdir) // '/tables/cx/cx_2_1_warmTarget.cdf'
       ier = NF90_OPEN(table_str,NF90_NOWRITE,ncid)
       ier = NF90_INQ_VARID(ncid,'axis_param_warmTarget',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(cx_2_1_axis(dimlen(1)))
       ier = NF90_GET_VAR(ncid, varid, cx_2_1_axis)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_INQ_VARID(ncid,'btsigv',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(cx_2_1_btsigv(dimlen(1),dimlen(2)))
       ier = NF90_GET_VAR(ncid, varid, cx_2_1_btsigv)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_CLOSE(ncid)

       ! Load Ion Impact CX 2-2
       table_str=TRIM(adasdir) // '/tables/cx/cx_2_2_warmTarget.cdf'
       ier = NF90_OPEN(table_str,NF90_NOWRITE,ncid)
       ier = NF90_INQ_VARID(ncid,'axis_param_warmTarget',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(cx_2_2_axis(dimlen(1)))
       ier = NF90_GET_VAR(ncid, varid, cx_2_2_axis)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_INQ_VARID(ncid,'btsigv',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(cx_2_2_btsigv(dimlen(1),dimlen(2)))
       ier = NF90_GET_VAR(ncid, varid, cx_2_2_btsigv)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_CLOSE(ncid)

       ! Load Ion Impact II 1-1
       table_str=TRIM(adasdir) // '/tables/ii/ii_1_1_warmTarget.cdf'
       ier = NF90_OPEN(table_str,NF90_NOWRITE,ncid)
       ier = NF90_INQ_VARID(ncid,'axis_param_warmTarget',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_1_1_axis(dimlen(1)))
       ier = NF90_GET_VAR(ncid, varid, ii_1_1_axis)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_INQ_VARID(ncid,'btsigv',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_1_1_btsigv(dimlen(1),dimlen(2)))
       ier = NF90_GET_VAR(ncid, varid, ii_1_1_btsigv)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_CLOSE(ncid)

       ! Load Ion Impact II 1-2
       table_str=TRIM(adasdir) // '/tables/ii/ii_1_2_warmTarget.cdf'
       ier = NF90_OPEN(table_str,NF90_NOWRITE,ncid)
       ier = NF90_INQ_VARID(ncid,'axis_param_warmTarget',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_1_2_axis(dimlen(1)))
       ier = NF90_GET_VAR(ncid, varid, ii_1_2_axis)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_INQ_VARID(ncid,'btsigv',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_1_2_btsigv(dimlen(1),dimlen(2)))
       ier = NF90_GET_VAR(ncid, varid, ii_1_2_btsigv)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_CLOSE(ncid)

       ! Load Ion Impact II 1-3
       table_str=TRIM(adasdir) // '/tables/ii/ii_1_3_warmTarget.cdf'
       ier = NF90_OPEN(table_str,NF90_NOWRITE,ncid)
       ier = NF90_INQ_VARID(ncid,'axis_param_warmTarget',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_1_3_axis(dimlen(1)))
       ier = NF90_GET_VAR(ncid, varid, ii_1_3_axis)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_INQ_VARID(ncid,'btsigv',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_1_3_btsigv(dimlen(1),dimlen(2)))
       ier = NF90_GET_VAR(ncid, varid, ii_1_3_btsigv)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_CLOSE(ncid)

       ! Load Ion Impact II 1-4
       table_str=TRIM(adasdir) // '/tables/ii/ii_1_4_warmTarget.cdf'
       ier = NF90_OPEN(table_str,NF90_NOWRITE,ncid)
       ier = NF90_INQ_VARID(ncid,'axis_param_warmTarget',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_1_4_axis(dimlen(1)))
       ier = NF90_GET_VAR(ncid, varid, ii_1_4_axis)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_INQ_VARID(ncid,'btsigv',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_1_4_btsigv(dimlen(1),dimlen(2)))
       ier = NF90_GET_VAR(ncid, varid, ii_1_4_btsigv)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_CLOSE(ncid)

       ! Load Ion Impact II 1-5
       table_str=TRIM(adasdir) // '/tables/ii/ii_1_5_warmTarget.cdf'
       ier = NF90_OPEN(table_str,NF90_NOWRITE,ncid)
       ier = NF90_INQ_VARID(ncid,'axis_param_warmTarget',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_1_5_axis(dimlen(1)))
       ier = NF90_GET_VAR(ncid, varid, ii_1_5_axis)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_INQ_VARID(ncid,'btsigv',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_1_5_btsigv(dimlen(1),dimlen(2)))
       ier = NF90_GET_VAR(ncid, varid, ii_1_5_btsigv)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_CLOSE(ncid)

       ! Load Ion Impact II 1-6
       table_str=TRIM(adasdir) // '/tables/ii/ii_1_6_warmTarget.cdf'
       ier = NF90_OPEN(table_str,NF90_NOWRITE,ncid)
       ier = NF90_INQ_VARID(ncid,'axis_param_warmTarget',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_1_6_axis(dimlen(1)))
       ier = NF90_GET_VAR(ncid, varid, ii_1_6_axis)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_INQ_VARID(ncid,'btsigv',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_1_6_btsigv(dimlen(1),dimlen(2)))
       ier = NF90_GET_VAR(ncid, varid, ii_1_6_btsigv)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_CLOSE(ncid)

       ! Load Ion Impact II 1-7
       table_str=TRIM(adasdir) // '/tables/ii/ii_1_7_warmTarget.cdf'
       ier = NF90_OPEN(table_str,NF90_NOWRITE,ncid)
       ier = NF90_INQ_VARID(ncid,'axis_param_warmTarget',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_1_7_axis(dimlen(1)))
       ier = NF90_GET_VAR(ncid, varid, ii_1_7_axis)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_INQ_VARID(ncid,'btsigv',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_1_7_btsigv(dimlen(1),dimlen(2)))
       ier = NF90_GET_VAR(ncid, varid, ii_1_7_btsigv)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_CLOSE(ncid)

       ! Load Ion Impact II 1-8
       table_str=TRIM(adasdir) // '/tables/ii/ii_1_8_warmTarget.cdf'
       ier = NF90_OPEN(table_str,NF90_NOWRITE,ncid)
       ier = NF90_INQ_VARID(ncid,'axis_param_warmTarget',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_1_8_axis(dimlen(1)))
       ier = NF90_GET_VAR(ncid, varid, ii_1_8_axis)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_INQ_VARID(ncid,'btsigv',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_1_8_btsigv(dimlen(1),dimlen(2)))
       ier = NF90_GET_VAR(ncid, varid, ii_1_8_btsigv)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_CLOSE(ncid)

       ! Load Ion Impact II 1-9
       table_str=TRIM(adasdir) // '/tables/ii/ii_1_9_warmTarget.cdf'
       ier = NF90_OPEN(table_str,NF90_NOWRITE,ncid)
       ier = NF90_INQ_VARID(ncid,'axis_param_warmTarget',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_1_9_axis(dimlen(1)))
       ier = NF90_GET_VAR(ncid, varid, ii_1_9_axis)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_INQ_VARID(ncid,'btsigv',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_1_9_btsigv(dimlen(1),dimlen(2)))
       ier = NF90_GET_VAR(ncid, varid, ii_1_9_btsigv)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_CLOSE(ncid)

       ! Load Ion Impact II 1-10
       table_str=TRIM(adasdir) // '/tables/ii/ii_1_10_warmTarget.cdf'
       ier = NF90_OPEN(table_str,NF90_NOWRITE,ncid)
       ier = NF90_INQ_VARID(ncid,'axis_param_warmTarget',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_1_10_axis(dimlen(1)))
       ier = NF90_GET_VAR(ncid, varid, ii_1_1_axis)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_INQ_VARID(ncid,'btsigv',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_1_10_btsigv(dimlen(1),dimlen(2)))
       ier = NF90_GET_VAR(ncid, varid, ii_1_10_btsigv)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_CLOSE(ncid)

       ! Load Ion Impact II 2-1
       table_str=TRIM(adasdir) // '/tables/ii/ii_2_1_warmTarget.cdf'
       ier = NF90_OPEN(table_str,NF90_NOWRITE,ncid)
       ier = NF90_INQ_VARID(ncid,'axis_param_warmTarget',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_2_1_axis(dimlen(1)))
       ier = NF90_GET_VAR(ncid, varid, ii_2_1_axis)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_INQ_VARID(ncid,'btsigv',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_2_1_btsigv(dimlen(1),dimlen(2)))
       ier = NF90_GET_VAR(ncid, varid, ii_2_1_btsigv)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_CLOSE(ncid)

       ! Load Ion Impact II 2-2
       table_str=TRIM(adasdir) // '/tables/ii/ii_2_2_warmTarget.cdf'
       ier = NF90_OPEN(table_str,NF90_NOWRITE,ncid)
       ier = NF90_INQ_VARID(ncid,'axis_param_warmTarget',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_2_2_axis(dimlen(1)))
       ier = NF90_GET_VAR(ncid, varid, ii_2_2_axis)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_INQ_VARID(ncid,'btsigv',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_2_2_btsigv(dimlen(1),dimlen(2)))
       ier = NF90_GET_VAR(ncid, varid, ii_2_2_btsigv)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_CLOSE(ncid)

       ! Load Ion Impact II 2-3
       table_str=TRIM(adasdir) // '/tables/ii/ii_2_3_warmTarget.cdf'
       ier = NF90_OPEN(table_str,NF90_NOWRITE,ncid)
       ier = NF90_INQ_VARID(ncid,'axis_param_warmTarget',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_2_3_axis(dimlen(1)))
       ier = NF90_GET_VAR(ncid, varid, ii_2_3_axis)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_INQ_VARID(ncid,'btsigv',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_2_3_btsigv(dimlen(1),dimlen(2)))
       ier = NF90_GET_VAR(ncid, varid, ii_2_3_btsigv)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_CLOSE(ncid)

       ! Load Ion Impact II 2-4
       table_str=TRIM(adasdir) // '/tables/ii/ii_2_4_warmTarget.cdf'
       ier = NF90_OPEN(table_str,NF90_NOWRITE,ncid)
       ier = NF90_INQ_VARID(ncid,'axis_param_warmTarget',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_2_4_axis(dimlen(1)))
       ier = NF90_GET_VAR(ncid, varid, ii_2_4_axis)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_INQ_VARID(ncid,'btsigv',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_2_4_btsigv(dimlen(1),dimlen(2)))
       ier = NF90_GET_VAR(ncid, varid, ii_2_4_btsigv)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_CLOSE(ncid)

       ! Load Ion Impact II 2-5
       table_str=TRIM(adasdir) // '/tables/ii/ii_2_5_warmTarget.cdf'
       ier = NF90_OPEN(table_str,NF90_NOWRITE,ncid)
       ier = NF90_INQ_VARID(ncid,'axis_param_warmTarget',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_2_5_axis(dimlen(1)))
       ier = NF90_GET_VAR(ncid, varid, ii_2_5_axis)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_INQ_VARID(ncid,'btsigv',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_2_5_btsigv(dimlen(1),dimlen(2)))
       ier = NF90_GET_VAR(ncid, varid, ii_2_5_btsigv)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_CLOSE(ncid)

       ! Load Ion Impact II 2-6
       table_str=TRIM(adasdir) // '/tables/ii/ii_2_6_warmTarget.cdf'
       ier = NF90_OPEN(table_str,NF90_NOWRITE,ncid)
       ier = NF90_INQ_VARID(ncid,'axis_param_warmTarget',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_2_6_axis(dimlen(1)))
       ier = NF90_GET_VAR(ncid, varid, ii_2_6_axis)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_INQ_VARID(ncid,'btsigv',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_2_6_btsigv(dimlen(1),dimlen(2)))
       ier = NF90_GET_VAR(ncid, varid, ii_2_6_btsigv)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_CLOSE(ncid)

       ! Load Ion Impact II 2-7
       table_str=TRIM(adasdir) // '/tables/ii/ii_2_7_warmTarget.cdf'
       ier = NF90_OPEN(table_str,NF90_NOWRITE,ncid)
       ier = NF90_INQ_VARID(ncid,'axis_param_warmTarget',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_2_7_axis(dimlen(1)))
       ier = NF90_GET_VAR(ncid, varid, ii_2_7_axis)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_INQ_VARID(ncid,'btsigv',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_2_7_btsigv(dimlen(1),dimlen(2)))
       ier = NF90_GET_VAR(ncid, varid, ii_2_7_btsigv)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_CLOSE(ncid)

       ! Load Ion Impact II 2-8
       table_str=TRIM(adasdir) // '/tables/ii/ii_2_8_warmTarget.cdf'
       ier = NF90_OPEN(table_str,NF90_NOWRITE,ncid)
       ier = NF90_INQ_VARID(ncid,'axis_param_warmTarget',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_2_8_axis(dimlen(1)))
       ier = NF90_GET_VAR(ncid, varid, ii_2_8_axis)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_INQ_VARID(ncid,'btsigv',varid)
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,NDIMS = ndims)
       ALLOCATE(dimid(ndims),dimlen(ndims))
       ier = NF90_INQUIRE_VARIABLE(ncid,varid,DIMIDS = dimid)
       DO i = 1, ndims
          ier = NF90_INQUIRE_DIMENSION(ncid,dimid(i),LEN = dimlen(i))
       END DO
       ALLOCATE(ii_2_8_btsigv(dimlen(1),dimlen(2)))
       ier = NF90_GET_VAR(ncid, varid, ii_2_8_btsigv)
       DEALLOCATE(dimid,dimlen)
       ier = NF90_CLOSE(ncid)
#endif

       END SUBROUTINE adas_load_tables

       SUBROUTINE deallocate_adas_tables
       IMPLICIT NONE
       IF (ALLOCATED(ei_1_axis)) DEALLOCATE(ei_1_axis)
       IF (ALLOCATED(ei_2_axis)) DEALLOCATE(ei_2_axis)
       IF (ALLOCATED(cx_1_1_axis)) DEALLOCATE(cx_1_1_axis)
       IF (ALLOCATED(cx_2_1_axis)) DEALLOCATE(cx_2_1_axis)
       IF (ALLOCATED(cx_2_2_axis)) DEALLOCATE(cx_2_2_axis)
       IF (ALLOCATED(ii_1_1_axis)) DEALLOCATE(ii_1_1_axis)
       IF (ALLOCATED(ii_1_2_axis)) DEALLOCATE(ii_1_2_axis)
       IF (ALLOCATED(ii_1_3_axis)) DEALLOCATE(ii_1_3_axis)
       IF (ALLOCATED(ii_1_4_axis)) DEALLOCATE(ii_1_4_axis)
       IF (ALLOCATED(ii_1_5_axis)) DEALLOCATE(ii_1_5_axis)
       IF (ALLOCATED(ii_1_6_axis)) DEALLOCATE(ii_1_6_axis)
       IF (ALLOCATED(ii_1_7_axis)) DEALLOCATE(ii_1_7_axis)
       IF (ALLOCATED(ii_1_8_axis)) DEALLOCATE(ii_1_8_axis)
       IF (ALLOCATED(ii_1_9_axis)) DEALLOCATE(ii_1_9_axis)
       IF (ALLOCATED(ii_1_10_axis)) DEALLOCATE(ii_1_10_axis)
       IF (ALLOCATED(ii_2_1_axis)) DEALLOCATE(ii_2_1_axis)
       IF (ALLOCATED(ii_2_2_axis)) DEALLOCATE(ii_2_2_axis)
       IF (ALLOCATED(ii_2_3_axis)) DEALLOCATE(ii_2_3_axis)
       IF (ALLOCATED(ii_2_4_axis)) DEALLOCATE(ii_2_4_axis)
       IF (ALLOCATED(ii_2_5_axis)) DEALLOCATE(ii_2_5_axis)
       IF (ALLOCATED(ii_2_6_axis)) DEALLOCATE(ii_2_6_axis)
       IF (ALLOCATED(ii_2_7_axis)) DEALLOCATE(ii_2_7_axis)
       IF (ALLOCATED(ii_2_8_axis)) DEALLOCATE(ii_2_8_axis)

       IF (ALLOCATED(ei_1_sigv)) DEALLOCATE(ei_1_sigv)
       IF (ALLOCATED(ei_2_sigv)) DEALLOCATE(ei_2_sigv)
       IF (ALLOCATED(cx_1_1_btsigv)) DEALLOCATE(cx_1_1_btsigv)
       IF (ALLOCATED(cx_2_1_btsigv)) DEALLOCATE(cx_2_1_btsigv)
       IF (ALLOCATED(cx_2_2_btsigv)) DEALLOCATE(cx_2_2_btsigv)
       IF (ALLOCATED(ii_1_1_btsigv)) DEALLOCATE(ii_1_1_btsigv)
       IF (ALLOCATED(ii_1_2_btsigv)) DEALLOCATE(ii_1_2_btsigv)
       IF (ALLOCATED(ii_1_3_btsigv)) DEALLOCATE(ii_1_3_btsigv)
       IF (ALLOCATED(ii_1_4_btsigv)) DEALLOCATE(ii_1_4_btsigv)
       IF (ALLOCATED(ii_1_5_btsigv)) DEALLOCATE(ii_1_5_btsigv)
       IF (ALLOCATED(ii_1_6_btsigv)) DEALLOCATE(ii_1_6_btsigv)
       IF (ALLOCATED(ii_1_7_btsigv)) DEALLOCATE(ii_1_7_btsigv)
       IF (ALLOCATED(ii_1_8_btsigv)) DEALLOCATE(ii_1_8_btsigv)
       IF (ALLOCATED(ii_1_9_btsigv)) DEALLOCATE(ii_1_9_btsigv)
       IF (ALLOCATED(ii_1_10_btsigv)) DEALLOCATE(ii_1_10_btsigv)
       IF (ALLOCATED(ii_2_1_btsigv)) DEALLOCATE(ii_2_1_btsigv)
       IF (ALLOCATED(ii_2_2_btsigv)) DEALLOCATE(ii_2_2_btsigv)
       IF (ALLOCATED(ii_2_3_btsigv)) DEALLOCATE(ii_2_3_btsigv)
       IF (ALLOCATED(ii_2_4_btsigv)) DEALLOCATE(ii_2_4_btsigv)
       IF (ALLOCATED(ii_2_5_btsigv)) DEALLOCATE(ii_2_5_btsigv)
       IF (ALLOCATED(ii_2_6_btsigv)) DEALLOCATE(ii_2_6_btsigv)
       IF (ALLOCATED(ii_2_7_btsigv)) DEALLOCATE(ii_2_7_btsigv)
       IF (ALLOCATED(ii_2_8_btsigv)) DEALLOCATE(ii_2_8_btsigv)
       END SUBROUTINE deallocate_adas_tables

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
       real*8  :: zion          !     izion        : atomic charge of secondary particle {1, ..., 10}
       integer :: npts_e, npts_t ! number of points for table
       real*8  :: xlr, xlr_t ! If XLR>0, then the X grid is equally spaced on a logarithmic scale:
                  ! LINEAR IF XLR .LE. 0.0
       real*8, dimension(n1) :: sigv_adas_wrk1,sigv_adas_wrk2 !sigma*v [m**3/sec]
       !checking input data for availability
       istat=0
       sigv_adas=0
       sigv_adas_wrk1 =0; sigv_adas_wrk2=0
       zion = DBLE(izneut_in)

       IF (freact_type == 1) THEN
          IF(beamchrg.eq.NEUTRAL) THEN
             xlr    =  log(cx_1_1_axis(2)/cx_1_1_axis(1))
             xlr_t  =  log(cx_1_1_axis(6)/cx_1_1_axis(5))
             npts_e = SIZE(cx_1_1_btsigv,DIM=1)
             npts_t = SIZE(cx_1_1_btsigv,DIM=2)
             call FLINT_Z(evec,tevec,sigv_adas_wrk1,n1,cx_1_1_btsigv, cx_1_1_axis(1),cx_1_1_axis(2),&
                                               xlr,npts_e,cx_1_1_axis(5),cx_1_1_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
             sigv_adas = sigv_adas_wrk1*(zion) ! This is the ADAS WAY
          ELSE IF (beamchrg.eq.ION) THEN
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
          IF(beamchrg.eq.NEUTRAL) THEN
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
          ELSE IF (beamchrg.eq.ION) THEN
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
             ELSE IF (zion >=2 .and. zion < 3) THEN
                xlr=log(ii_2_2_axis(2)/ii_2_2_axis(1))
                xlr_t=log(ii_2_2_axis(6)/ii_2_2_axis(5))
                npts_e = SIZE(ii_2_2_btsigv,DIM=1)
                npts_t = SIZE(ii_2_2_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk1,n1,ii_2_2_btsigv,ii_2_2_axis(1),ii_2_2_axis(2),&
                             xlr,npts_e,ii_2_2_axis(5),ii_2_2_axis(6),&
                             xlr_t,npts_t,iwarn_adas)
                xlr=log(ii_2_3_axis(2)/ii_2_3_axis(1))
                xlr_t=log(ii_2_3_axis(6)/ii_2_3_axis(5))
                npts_e = SIZE(ii_2_3_btsigv,DIM=1)
                npts_t = SIZE(ii_2_3_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk2,n1,ii_2_3_btsigv,ii_2_3_axis(1),ii_2_3_axis(2),&
                             xlr,npts_e,ii_2_3_axis(5),ii_2_3_axis(6),&
                             xlr_t,npts_t,iwarn_adas)
                sigv_adas = sigv_adas_wrk1 + (sigv_adas_wrk2-sigv_adas_wrk1)*(zion-2)
             ELSE IF (zion >=3 .and. zion < 4) THEN
                xlr    =  log(ii_2_3_axis(2)/ii_2_3_axis(1))
                xlr_t  =  log(ii_2_3_axis(6)/ii_2_3_axis(5))
                npts_e = SIZE(ii_2_3_btsigv,DIM=1)
                npts_t = SIZE(ii_2_3_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk1,n1,ii_2_3_btsigv, ii_2_3_axis(1),ii_2_3_axis(2),&
                                               xlr,npts_e,ii_2_3_axis(5),ii_2_3_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
                xlr    =  log(ii_2_4_axis(2)/ii_2_4_axis(1))
                xlr_t  =  log(ii_2_4_axis(6)/ii_2_4_axis(5))
                npts_e = SIZE(ii_2_4_btsigv,DIM=1)
                npts_t = SIZE(ii_2_4_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk2,n1,ii_2_4_btsigv, ii_2_4_axis(1),ii_2_4_axis(2),&
                                               xlr,npts_e,ii_2_4_axis(5),ii_2_4_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
                sigv_adas = sigv_adas_wrk1 + (sigv_adas_wrk2-sigv_adas_wrk1)*(zion-3)
             ELSE IF (zion >=4 .and. zion < 5) THEN
                xlr    =  log(ii_2_4_axis(2)/ii_2_4_axis(1))
                xlr_t  =  log(ii_2_4_axis(6)/ii_2_4_axis(5))
                npts_e = SIZE(ii_2_4_btsigv,DIM=1)
                npts_t = SIZE(ii_2_4_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk1,n1,ii_2_4_btsigv, ii_2_4_axis(1),ii_2_4_axis(2),&
                                               xlr,npts_e,ii_2_4_axis(5),ii_2_4_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
                xlr    =  log(ii_2_5_axis(2)/ii_2_5_axis(1))
                xlr_t  =  log(ii_2_5_axis(6)/ii_2_5_axis(5))
                npts_e = SIZE(ii_2_5_btsigv,DIM=1)
                npts_t = SIZE(ii_2_5_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk2,n1,ii_2_5_btsigv, ii_2_5_axis(1),ii_2_5_axis(2),&
                                               xlr,npts_e,ii_2_5_axis(5),ii_2_5_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
                sigv_adas = sigv_adas_wrk1 + (sigv_adas_wrk2-sigv_adas_wrk1)*(zion-4)
             ELSE IF (zion >=5 .and. zion < 6) THEN
                xlr    =  log(ii_2_5_axis(2)/ii_2_5_axis(1))
                xlr_t  =  log(ii_2_5_axis(6)/ii_2_5_axis(5))
                npts_e = SIZE(ii_2_5_btsigv,DIM=1)
                npts_t = SIZE(ii_2_5_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk1,n1,ii_2_5_btsigv, ii_2_5_axis(1),ii_2_5_axis(2),&
                                               xlr,npts_e,ii_2_5_axis(5),ii_2_5_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
                xlr    =  log(ii_2_6_axis(2)/ii_2_6_axis(1))
                xlr_t  =  log(ii_2_6_axis(6)/ii_2_6_axis(5))
                npts_e = SIZE(ii_2_6_btsigv,DIM=1)
                npts_t = SIZE(ii_2_6_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk2,n1,ii_2_6_btsigv, ii_2_6_axis(1),ii_2_6_axis(2),&
                                               xlr,npts_e,ii_2_6_axis(5),ii_2_6_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
                sigv_adas = sigv_adas_wrk1 + (sigv_adas_wrk2-sigv_adas_wrk1)*(zion-5)
             ELSE IF (zion >=6 .and. zion < 7) THEN
                xlr    =  log(ii_2_6_axis(2)/ii_2_6_axis(1))
                xlr_t  =  log(ii_2_6_axis(6)/ii_2_6_axis(5))
                npts_e = SIZE(ii_2_6_btsigv,DIM=1)
                npts_t = SIZE(ii_2_6_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk1,n1,ii_2_6_btsigv, ii_2_6_axis(1),ii_2_6_axis(2),&
                                               xlr,npts_e,ii_2_6_axis(5),ii_2_6_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
                xlr    =  log(ii_2_7_axis(2)/ii_2_7_axis(1))
                xlr_t  =  log(ii_2_7_axis(6)/ii_2_7_axis(5))
                npts_e = SIZE(ii_2_7_btsigv,DIM=1)
                npts_t = SIZE(ii_2_7_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk2,n1,ii_2_7_btsigv, ii_2_7_axis(1),ii_2_7_axis(2),&
                                               xlr,npts_e,ii_2_7_axis(5),ii_2_7_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
                sigv_adas = sigv_adas_wrk1 + (sigv_adas_wrk2-sigv_adas_wrk1)*(zion-6)
             ELSE IF (zion >=7 .and. zion < 8) THEN
                xlr    =  log(ii_2_7_axis(2)/ii_2_7_axis(1))
                xlr_t  =  log(ii_2_7_axis(6)/ii_2_7_axis(5))
                npts_e = SIZE(ii_2_7_btsigv,DIM=1)
                npts_t = SIZE(ii_2_7_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk1,n1,ii_2_7_btsigv, ii_2_7_axis(1),ii_2_7_axis(2),&
                                               xlr,npts_e,ii_2_7_axis(5),ii_2_7_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
                xlr    =  log(ii_2_8_axis(2)/ii_2_8_axis(1))
                xlr_t  =  log(ii_2_8_axis(6)/ii_2_8_axis(5))
                npts_e = SIZE(ii_2_8_btsigv,DIM=1)
                npts_t = SIZE(ii_2_8_btsigv,DIM=2)
                call FLINT_Z(evec,tevec,sigv_adas_wrk2,n1,ii_2_8_btsigv, ii_2_8_axis(1),ii_2_8_axis(2),&
                                               xlr,npts_e,ii_2_8_axis(5),ii_2_8_axis(6),&
                                               xlr_t,npts_t,iwarn_adas)
                sigv_adas = sigv_adas_wrk1 + (sigv_adas_wrk2-sigv_adas_wrk1)*(zion-7)
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

end module adas_mod_parallel