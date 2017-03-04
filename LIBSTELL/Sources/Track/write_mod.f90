      MODULE WRITE_MOD
!-------------------------------------------------------------------------------
!WRITE-WRITEs output in standardized formats
!
!WRITE_MOD is an F90 module of standarized output routines
!
!References:
!
! W.A.Houlberg 7/2001
!
!Contains PUBLIC routines:
!
!  WRITE_OUT1          -writes 1D data
!  WRITE_C             -writes character variables
!  WRITE_IR            -writes integer and real variables
!  WRITE_LINE          -writes a line
!  WRITE_LINE_IR       -writes a writes a line plus integer and real numbers
!
!Comments:
!
!  The modernization of the code structure into an F90 module takes advantage of
!    some of the more attractive features of F90:
!    -use of KIND for precision declarations
!    -optional arguments for I/O
!    -generic names for all intrinsic functions
!    -compilation using either free or fixed form
!    -no common blocks or other deprecated Fortran features
!    -dynamic and automatic alocation of variables
!    -array syntax for vector operations
!-------------------------------------------------------------------------------
      USE SPEC_KIND_MOD
      IMPLICIT NONE

!-------------------------------------------------------------------------------
! Public procedures
!-------------------------------------------------------------------------------
      CONTAINS

      SUBROUTINE WRITE_OUT1(nout,c_nout,nx,nv,values,names,units,k_lr,         &
     &                      LABEL)
!-------------------------------------------------------------------------------
!WRITE_OUT1 writes 1D data to one of several files
!References:
!  W.A.Houlberg 7/2001
!Input:
!  nout                -output file unit number [-]
!  c_nout              -type of output file [character]
!                      -'1d' for graphical post-processing file
!                      -'sum' for summary file
!                      -'netcdf' for netcdf file
!  nx                  -no. of variables [-]
!  nv                  -no. of points [-]
!  values(nx,nv)       -values [-]
!  names(nv)           -names of variables [character]
!  units(nv)           -units of variables [character]
!  k_lr                -option for shifting names and units [-]
!                      =-1 shift left (remove leading blanks)
!                      =+1 shift right (insert leading blanks)
!                      =else no shift
!Optional input:
!  LABEL               -label to insert as break in output file
!-------------------------------------------------------------------------------

!Declaration of input variables
      CHARACTER(len=*), INTENT(IN) ::                                          &
     & c_nout,names(:),units(:)

      INTEGER, INTENT(IN) ::                                                   &
     & k_lr,nout,nx,nv

      REAL(KIND=rspec), INTENT(IN) ::                                          &
     & values(:,:)

!Declaration of optional input variables
      CHARACTER(len=*), INTENT(IN), OPTIONAL ::                                &
     & LABEL

!Declaration of local variables
      CHARACTER(len=120) ::                                                    &
     & l

      CHARACTER(len=15) ::                                                     &
     & cdum(6)

      CHARACTER(len=15), ALLOCATABLE ::                                        &
     & n(:),u(:)

      INTEGER ::                                                               &
     & i,j,nx5,nv5,idum(5)

      REAL(KIND=rspec) ::                                                      &
     & rdum(5)

      REAL(KIND=rspec), ALLOCATABLE ::                                         &
     & v(:,:)

!-------------------------------------------------------------------------------
!Initialization
!-------------------------------------------------------------------------------
!Pad arrays to an even multiple of 5
      nx5=((nx+4)/5)*5
      nv5=((nv+4)/5)*5

!Allocate padded temporary arrays
      ALLOCATE(n(1:nv5),u(1:nv5),v(1:nx5,1:nv5))
      n(:)='dummy'
      u(:)='-'
      v(:,:)=0.0_rspec

!Copy and pad info
      n(1:nv)=names(1:nv)
      u(1:nv)=units(1:nv)
      v(1:nx,1:nv)=values(1:nx,1:nv)

!Shift names and units right, left, or leave untouched

      IF(k_lr == -1) THEN

        n(:)=ADJUSTL(n)
        u(:)=ADJUSTL(u)

      ELSEIF(k_lr == 1) THEN

        n(:)=ADJUSTR(n)
        u(:)=ADJUSTR(u)

      ENDIF
      
!-------------------------------------------------------------------------------
!Print data to 1D file
!-------------------------------------------------------------------------------
      IF(c_nout == '1d') THEN

!Number of profiles and radial nodes
        idum(1)=nv
        idum(2)=nx
        rdum(1)=0.0_rspec
        CALL WRITE_IR(nout,2,idum,0,rdum,15,0)

!Names
        DO j=1,nv5/5
 
          CALL WRITE_C(nout,5,n(5*j-4),15)

        ENDDO

!Units
        DO j=1,nv5/5

          CALL WRITE_C(nout,5,u(5*j-4),15)

        ENDDO

!Values
        idum(:)=0

        DO j=1,nv5

          DO i=1,nx5/5

            CALL WRITE_IR(nout,0,idum,5,v(1+5*(i-1),j),15,2)

          ENDDO

        ENDDO

      ENDIF
      
!-------------------------------------------------------------------------------
!Print data to summary file
!-------------------------------------------------------------------------------
      IF(c_nout == 'sum') THEN

!Label
        IF(PRESENT(LABEL)) CALL WRITE_LINE(nout,LABEL,1,1)

        l=''

        DO j=1,nv5/5

!Names
          cdum(1)='              i'
          cdum(2:6)=n(5*(j-1)+1:5)
          CALL WRITE_LINE(nout,l,0,0)
          CALL WRITE_C(nout,6,cdum,15)

!Units
          cdum(1)='              -'
          cdum(2:6)=u(5*(j-1)+1:5)
          CALL WRITE_C(nout,6,cdum,15)

!Values
          DO i=1,nx

            idum(1)=i
            rdum(1:5)=v(i,5*(j-1)+1:5)
            CALL WRITE_IR(nout,1,idum,5,rdum,15,2)

          ENDDO

        ENDDO

      ENDIF

      END SUBROUTINE WRITE_OUT1

      SUBROUTINE WRITE_C(nout,n_c,c,n_l)
!-------------------------------------------------------------------------------
!WRITE_C writes out n_c character variables in fields of length n_l
!References:
!  W.A.Houlberg 7/2001
!Input:
!  nout                -output file unit number [-]
!  n_c                 -number of character variables [-]
!  c(n_c)              -column headings [character]
!  n_l                 -length of field [-]
!Comments:
!  This routine can be used to write out a set of column headings, variable
!    names or other appications that use a set of character strings
!-------------------------------------------------------------------------------

!Declaration of input variables
      INTEGER, INTENT(IN) ::                                                   &
     & n_c,                                                                    &
     & n_l,                                                                    &
     & nout

      CHARACTER(len=*), INTENT(IN) ::                                          &
     & c(1:n_c)

!Declaration of local variables
      CHARACTER(len=30) ::                                                     &
     & char

      CHARACTER(len=2) ::                                                      &
     & c_c,c_l

      INTEGER ::                                                               &
     & i

!Use internal write to set number and length of fields
      WRITE(c_c,'(i2)') n_c
      WRITE(c_l,'(i2)') n_l

!Set format
      char='(1x,'//c_c//'a'//c_l//')'

!Write 
      WRITE(nout,char) (c(i),i=1,n_c)

      END SUBROUTINE WRITE_C

      SUBROUTINE WRITE_IR(nout,n_i,i,n_r,r,n_l,k_format)
!-------------------------------------------------------------------------------
!WRITE_IR writes n_i integer variables followed by n_r real variables
!References:
!  W.A.Houlberg 7/2001
!Input:
!  nout                -output file unit number [-]
!  n_i                 -number of integer variables [-]
!  i(n_i)              -integer array [-]
!  n_r                 -number of real variables [-]
!  r(n_r)              -real array [arb]
!  n_l                 -length of field [-]
!  k_format            -real format option [-]
!                      =1 use f
!                      =2 use 1pe
!                      =else use e
!-------------------------------------------------------------------------------

!Declaration of input variables
      INTEGER, INTENT(IN) ::                                                   &
     & k_format,                                                               &
     & n_i,                                                                    &
     & n_l,                                                                    &
     & n_r,                                                                    &
     & nout,                                                                   &
     & i(1:n_i)

      REAL(KIND=rspec), INTENT(IN) ::                                          &
     & r(1:n_r)

!Declaration of local variables
      CHARACTER(len=30) ::                                                     &
     & char

      CHARACTER(len=2)                                                         &
     & c_i,c_l,cp6_l,cp8_l,c_r

      INTEGER ::                                                               &
     & j

!Use internal write to set number and length of fields
      WRITE(c_i,'(i2)') n_i
      WRITE(c_l,'(i2)') n_l
      IF(c_l(1:1) == ' ') c_l(1:1)='0'
      WRITE(cp6_l,'(i2)') n_l-6
      IF(cp6_l(1:1) == ' ') cp6_l(1:1)='0'
      WRITE(cp8_l,'(i2)') n_l-8
      IF(cp8_l(1:1) == ' ') cp8_l(1:1)='0'
      WRITE(c_r,'(i2)') n_r

!Set format
      IF(n_i == 0) THEN

        !Real data only
        IF(k_format == 1) THEN

          !f format
          char='(1x,'//c_r//'(f'//c_l//'.'//cp6_l//'))'

        ELSEIF(k_format == 2) THEN

          !1pe format
          char='(1x,'//c_r//'(1pe'//c_l//'.'//cp8_l//'))'

        ELSE

          !e format
          char='(1x,'//c_r//'(e'//c_l//'.'//cp8_l//'))'

        ENDIF

        WRITE(nout,char) (r(j),j=1,n_r)

      ELSEIF(n_r == 0) THEN

        !Integer data only
        char='(1x,'//c_i//'(i'//c_l//'))'
        WRITE(nout,char) (i(j),j=1,n_i)

      ELSE

        !Both integer and real data
        IF(k_format == 1) THEN

          !f format
          char='(1x,'//c_i//'(i'//c_l//'),'//c_r                               &
     &         //'(f'//c_l//'.'//cp6_l//'))'

        ELSEIF(k_format == 2) THEN

          !1pe format
          char='(1x,'//c_i//'(i'//c_l//'),'//c_r                               &
     &         //'(1pe'//c_l//'.'//cp8_l//'))'

        ELSE

          !e format
          char='(1x,'//c_i//'(i'//c_l//'),'//c_r                               &
     &         //'(e'//c_l//'.'//cp8_l//'))'

        ENDIF

!Output line
        WRITE(nout,char) (i(j),j=1,n_i),(r(j),j=1,n_r)

      ENDIF

      END SUBROUTINE WRITE_IR

      SUBROUTINE WRITE_LINE(nout,label,k_above,k_below)
!-------------------------------------------------------------------------------
!WRITE_LINE writes a line (character string) preceeded by k_above blank lines
!  and followed by k_below blank lines
!References:
!  W.A.Houlberg 7/2001
!Input:
!  nout                -output file unit number [-]
!  label               -label to be printed [character]
!  k_above             -number of blank lines above label [-]
!  k_below             -number of blanklines below label [-]
!-------------------------------------------------------------------------------

!Declaration of input variables
      CHARACTER(len=*), INTENT(IN) ::                                          &
     & label

      INTEGER, INTENT(IN) ::                                                   &
     & k_above,k_below,                                                        &
     & nout

!Declaration of local variables
      INTEGER ::                                                               &
     & j

!Blank lines before text
      IF(k_above > 0) THEN

        DO j=1,k_above !Over leading lines

          WRITE(nout,'( )')

        ENDDO !Over leading lines

      ENDIF

!Text line
      WRITE(nout,'(a)') label

!Blank lines after text
      IF(k_below > 0) THEN

        DO j=1,k_below !Over trailing lines

          WRITE(nout,'( )')

        ENDDO !Over trailing lines

      ENDIF

      END SUBROUTINE WRITE_LINE

      SUBROUTINE WRITE_LINE_IR(nout,label,n_i,i,n_r,r,n_l,k_format)
!-------------------------------------------------------------------------------
!WRITE_LINE_IR writes a line followed by n_i integer numbers and n_r real
!  numbers
!References:
!  W.A.Houlberg 7/2001
!Input:
!  nout                -output file unit number [-]
!  label               -label to be printed [character]
!  n_i                 -number of integer variables [-]
!  i(n_i)              -integer array [-]
!  n_r                 -number of real variables [-]
!  r(n_r)              -real array [-]
!  n_l                 -length of field [-]
!  k_format            -format option for real variables [-]
!                      =1 use f
!                      =2 use 1pe
!                      =else use e
!-------------------------------------------------------------------------------

!Declaration of input variables
      CHARACTER(len=*), INTENT(IN) ::                                          &
     & label

      INTEGER, INTENT(IN) ::                                                   &
     & k_format,                                                               &
     & n_i,n_l,n_r,                                                            &
     & nout,                                                                   &
     & i(1:n_i)

      REAL(KIND=rspec), INTENT(IN) ::                                          &
     & r(1:n_r)

!Declaration of local variables
      CHARACTER(len=30) ::                                                     &
     & char

      CHARACTER(len=2) ::                                                      &
     & c_i,c_l,c_r

      INTEGER ::                                                               &
     & j

!Use internal write to set number and length of fields
      WRITE(c_i,'(i2)') n_i
      WRITE(c_r,'(i2)') n_r
      WRITE(c_l,'(i2)') n_l

!Set format
      IF(n_i == 0) THEN

        !Real data only
        IF(k_format == 1) THEN

          !f format
          char='(a48,'//c_r//'(f'//c_l//'.6))'

        ELSEIF(k_format == 2) THEN

          !1pe format
          char='(a48,'//c_r//'(1pe'//c_l//'.4))'

        ELSE

          !e format
          char='(a48,'//c_r//'(e'//c_l//'.4))'

        ENDIF

        WRITE(nout,char) label,(r(j),j=1,n_r)

      ELSEIF(n_r == 0) THEN

        !Integer data only
        char='(a48,'//c_i//'(i'//c_l//'))'
        WRITE(nout,char) label,(i(j),j=1,n_i)

      ELSE

        !Both integer and real data
        IF(k_format  == 1) THEN

          !f format
          char='(a48,'//c_i//'(i12),'//c_r//'(f'//c_l//'.6))'

        ELSEIF(k_format == 2) THEN

          !1pe format
          char='(a48,'//c_i//'(i12),'//c_r//'(1pe'//c_l//'.4))'

        ELSE

          !e format
          char='(a48,'//c_i//'(i12),'//c_r//'(e'//c_l//'.4))'

        ENDIF

!Output line
        WRITE(nout,char) label,(i(j),j=1,n_i),(r(j),j=1,n_r)

      ENDIF

      END SUBROUTINE WRITE_LINE_IR

      END MODULE WRITE_MOD
