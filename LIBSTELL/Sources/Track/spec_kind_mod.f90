MODULE SPEC_KIND_MOD
!-------------------------------------------------------------------------------
!For REAL variables:
!  Use p=6, r=35   for single precision on 32 bit machine
!  Use p=12,r=100  for double precision on 32 bit machine
!                  for single precision on 64 bit machine
!Parameters for SELECTED_REAL_KIND:
!  p                   -number of digits
!  r                   -range of exponent from -r to r
!-------------------------------------------------------------------------------
INTEGER, PARAMETER :: &
  rspec = SELECTED_REAL_KIND(p=12,r=100), &
  ispec = SELECTED_INT_KIND(r=8)

END MODULE SPEC_KIND_MOD
