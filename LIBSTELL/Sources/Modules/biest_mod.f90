MODULE biest_mod

USE iso_c_binding

IMPLICIT NONE

INTERFACE

SUBROUTINE biest_init(fname) BIND(C,name="biest_init")
  IMPORT
  CHARACTER(c_char) :: fname(*)
END SUBROUTINE

SUBROUTINE biest_finalize() BIND(C,name="biest_finalize")
  IMPORT
END SUBROUTINE

SUBROUTINE biest_fill_grid( &
    nr,nphi,nz, &
    R,PHI,Z, &
    BR,BPHI,BZ ) &
    BIND(C,name="biest_fill_grid")

  IMPORT

  INTEGER(c_int), VALUE :: nr,nphi,nz

  REAL(c_double) :: R(*)
  REAL(c_double) :: PHI(*)
  REAL(c_double) :: Z(*)

  REAL(c_double) :: BR(*)
  REAL(c_double) :: BPHI(*)
  REAL(c_double) :: BZ(*)

END SUBROUTINE

END INTERFACE

END MODULE
