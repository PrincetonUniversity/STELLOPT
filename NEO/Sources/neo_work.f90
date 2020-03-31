
MODULE neo_work
! Working parameters
  USE neo_precision

  REAL(kind=dp) ::   theta_start
  REAL(kind=dp) ::   theta_end
  REAL(kind=dp) ::   theta_int

  REAL(kind=dp) ::   phi_start
  REAL(kind=dp) ::   phi_end
  REAL(kind=dp) ::   phi_int

  REAL(kind=dp),    DIMENSION(:,:),     ALLOCATABLE :: cosmth,sinmth
  REAL(kind=dp),    DIMENSION(:,:),     ALLOCATABLE :: cosnph,sinnph
  REAL(kind=dp),    DIMENSION(:),       ALLOCATABLE :: theta_arr,phi_arr
  REAL(kind=dp),    DIMENSION(:,:),     ALLOCATABLE :: r,z,l,b
  REAL(kind=dp),    DIMENSION(:,:),     ALLOCATABLE :: r_tb,z_tb,p_tb,b_tb
  REAL(kind=dp),    DIMENSION(:,:),     ALLOCATABLE :: r_pb,z_pb,p_pb,b_pb
  REAL(kind=dp),    DIMENSION(:,:),     ALLOCATABLE :: gtbtb,gpbpb,gtbpb
  REAL(kind=dp),    DIMENSION(:,:),     ALLOCATABLE :: isqrg,sqrg11,kg,pard
  REAL(kind=dp),    DIMENSION(:,:),     ALLOCATABLE :: bqtphi
END MODULE neo_work
