!-----------------------------------------------------------------------
!     Subroutine:    gaunt_freefree
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          04/18/2021
!     Description:   This subroutine returns the Gaunt Free_Free Factor
!                    <g_{ff}> based on the paper
!                    Sutherland, R. S. (1998). Accurate free-free Gaunt
!                         factors for astrophysical plasmas. Monthly 
!                         Notices of the Royal Astronomical Society, 
!                         300(2), 321–330. 
!                        http://doi.org/10.1046/j.1365-8711.1998.01687.x
!                    It is based on interpolation provided for Table 2.
!                    It takes arguments
!                        gamma2: gamma^2 = Z*Z*Ry/(k*T)
!                        u:            u = hc/(lambda*k*T)
!-----------------------------------------------------------------------
      SUBROUTINE gaunt_freefree(gamma2,u,gff)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
!-----------------------------------------------------------------------
!     Arguments
!          gamma2       Gamma^2=Z*Z*Ry/(k*T)
!          u            u=hc/lambda*k*T
!          gff          Interpolated Gaunt Free-Free Factor
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), INTENT(in) :: gamma2
      REAL(rprec), INTENT(in) :: u
      REAL(rprec), INTENT(out) :: gff
!-----------------------------------------------------------------------
!     PARAMETERS
!          A     Table 2
!          A1    Rows of Table 2
!          A2    Columns of Table 2
!          Note the notation for the table is reversed.
!-----------------------------------------------------------------------
      REAL(rprec), DIMENSION(17,9) :: A
      REAL(rprec), DIMENSION(17), PARAMETER :: A1 = &
        (/-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0/)
      REAL(rprec), DIMENSION(9), PARAMETER :: A2 = &
        (/-4.0,-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0/)
!-----------------------------------------------------------------------
!     Variables
!          log10g2, log10u     Log base 10 of input quantities.
!          idex1,idex2         Index for spline tables
!          idex1b,idex2b       idex + 1
!          delta1,delta2       position between table units from [0,1]
!-----------------------------------------------------------------------
      INTEGER :: idex1,idex2, idex1b, idex2b
      REAL(rprec) :: log10g2, log10u, delta1, delta2
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      gff=0.0
      ! Intializations
      xlog10 = log10(x)
      idex = MIN(MAX(COUNT(table_log_g2<xlog10),1),41)
      delta = xlog10 - table_gff0(idex)
      gff = ( ( delta * table_s3(idex) + table_s2(idex) ) * delta + &
                                         table_s1(idex) ) * delta + &
            table_gff0(idex)
      RETURN

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE gaunt_freefree