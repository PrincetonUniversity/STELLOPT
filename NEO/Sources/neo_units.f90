
MODULE neo_units
! Units and Formats
  USE neo_precision
  INTEGER, PARAMETER ::   r_u2  = 4
  INTEGER, PARAMETER ::   r_us  = 5
  INTEGER, PARAMETER ::   r_u23 = 23
  INTEGER, PARAMETER ::   r_ua  = 21
  INTEGER, PARAMETER ::   w_u4  = 10
  INTEGER, PARAMETER ::   r_u10 = 3
  INTEGER, PARAMETER ::   w_us  = 6
  INTEGER, PARAMETER ::   w_u10 = 7
  INTEGER, PARAMETER ::   w_u20 = 8
  INTEGER, PARAMETER ::   w_u30 = 9
  INTEGER, PARAMETER ::   w_u50 = 11
  INTEGER, PARAMETER ::   w_u60 = 12
  INTEGER, PARAMETER ::   w_u70 = 13
  INTEGER, PARAMETER ::   w_u80 = 14
  INTEGER, PARAMETER ::   w_u90 = 15

  INTEGER            ::   w_u6_open, w_u1 = w_u10, w_u2 = w_u20, w_u3 = w_u30,      &
                          w_u5 = w_u50,w_u6 = w_u60, w_u7 = w_u70, w_u8 = w_u80,    &
                          w_u9 = w_u90, r_u1 = r_u10

  CHARACTER(20),PARAMETER :: format220="(500d18.5)"

  CHARACTER*(120) :: arg1, extension    ! LPK
  INTEGER  :: numargs                   ! LPK

END MODULE neo_units
